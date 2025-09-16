
#ifndef STREAMING_KMER_HASHSET_HPP
#define STREAMING_KMER_HASHSET_HPP



#include "DNA_Utility.hpp"
#include "Kmer.hpp"
#include "Directed_Vertex.hpp"
#include "Minimizer_Iterator.hpp"
#include "Hash_Table.hpp"
#include "Parser.hpp"
#include "utility.hpp"

#include "unordered_dense/unordered_dense.h"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <cassert>


// Hashset particularly suited for streaming k-mer operations.
template <uint16_t k, uint16_t l = 19>
class Streaming_Kmer_Hashset
{
private:

    static constexpr double lf_default = 0.8;   // Default load factor: 80%.
    static constexpr std::size_t B = 32;    // Bucket size.

    std::size_t cap_;   // Capacity of the table in terms of buckets; adjusted to be a power of 2.

    Kmer<k>* K; // Flat table of key-collection of size `cap x B`.
    uint8_t* S; // Size of each bucket.

    static constexpr uint64_t skew_th = B;
    Hash_Set<uint64_t> skew_mins;   // Skew minimizers.
    Hash_Set<Kmer<k>> overflow; // The overflow table.

    uint64_t sum_len;   // Total sequence length provided.
    uint64_t sup_kmer_c;    // Number of super k-mers in the provided sequences to the set.

    // Following are for used for developer stats.
    ankerl::unordered_dense::map<uint64_t, ankerl::unordered_dense::set<Kmer<k>, Kmer_Hasher<k>>> min_to_kmer;  // Set of overflowing minimizers.


    // Tries to insert the k-mer `key` into the bucket `b`, returning `true` iff
    // `key` is absent at bucket `b` after this operation.
    bool insert_at(Kmer<k> key, std::size_t b);

    // Returns `true` iff the k-mer `key` is present at the bucket `b`.
    bool find(Kmer<k> key, std::size_t b) const;

public:

    // Constructs a streaming k-mer hashset to support up-to `max_sz` k-mers,
    // with a maximum load-factor of `lf`.
    Streaming_Kmer_Hashset(std::size_t max_sz, double lf = lf_default);

    Streaming_Kmer_Hashset(const Streaming_Kmer_Hashset&) = delete;
    Streaming_Kmer_Hashset(Streaming_Kmer_Hashset&&) = delete;
    Streaming_Kmer_Hashset& operator=(const Streaming_Kmer_Hashset&) = delete;
    Streaming_Kmer_Hashset& operator=(Streaming_Kmer_Hashset&&) = delete;

    ~Streaming_Kmer_Hashset();

    // Returns the capacity (in terms of keys) of the table.
    auto capacity() const { return cap_ * B; }

    // Returns the size of the table.
    auto size() const { uint64_t sz = 0; std::for_each(S, S + cap_, [&](auto s){ sz += s; }); return sz + overflow.size(); }

    // Clears the table.
    void clear();

    // Tries to insert all the k-mers from the `l`-length sequence `S` into the
    // table.
    template <bool track_stats_ = false> void insert(const char* s, std::size_t len);

    // Preprocesses a prefix of sequences from the files at `file_path` to
    // populate a collection of skew-minimizers.
    void preprocess_skew_mins(const std::vector<std::string>& file_paths);

    // Reports developer statistics about the table.
    void report_stats() const;
};


template <uint16_t k, uint16_t l>
inline Streaming_Kmer_Hashset<k, l>::Streaming_Kmer_Hashset(const std::size_t max_sz, const double lf):
      cap_(ceil_pow_2((max_sz / lf) / B))
    , K(allocate<Kmer<k>>(cap_ * B))
    , S(allocate<uint8_t>(cap_))
{
    clear();
}


template <uint16_t k, uint16_t l>
inline Streaming_Kmer_Hashset<k, l>::~Streaming_Kmer_Hashset()
{
    deallocate(K), deallocate(S);
}


template <uint16_t k, uint16_t l>
inline void Streaming_Kmer_Hashset<k, l>::clear()
{
    sum_len = 0, sup_kmer_c = 0;

    std::memset(S, 0, cap_);
    skew_mins.clear();
    overflow.clear();

    min_to_kmer.clear();
}


template <uint16_t k, uint16_t l>
template <bool track_stats_>
inline void Streaming_Kmer_Hashset<k, l>::insert(const char* const s, const std::size_t len)
{
    sum_len += len;
    if(len < k)
        return;

    Directed_Vertex<k> v(s);
    std::size_t kmer_idx = 0;

    Min_Iterator<k, l> min_it;
    min_it.reset(s);
    uint64_t last_min = 0;  // Undercounting the first super k-mer from `s` with probability 1 / 2^64.
    while(true)
    {
        const auto key = v.canonical();
        const auto h = min_it.hash();
        const auto idx_mask = cap_ - 1;

        if constexpr(track_stats_)
        {
            sup_kmer_c += (last_min != h);
            last_min = h;
            min_to_kmer[h].insert(key);
        }

        if(skew_mins.contains(h))
            overflow.insert(key);
        else
        {
            const auto b = h & idx_mask;
            if(!insert_at(key, b))
            {
                const auto h_1 = XXH3_64bits_withSeed(&h, sizeof(h), 0);
                const auto h_2 = XXH3_64bits_withSeed(&h, sizeof(h), 1llu << 63);

                const auto p = h_1 & idx_mask;
                const auto q = h_2 & idx_mask;
                if(!find(key, p) && !find(key, q))
                {
                    if(S[p] < S[q])
                        K[p * B + S[p]++] = key;
                    else if(S[q] < B)
                        K[q * B + S[q]++] = key;
                    else
                        overflow.insert(key);
                }
            }
        }

        if(kmer_idx + k == len)
            break;

        min_it.advance(s[kmer_idx + k]);
        v.roll_forward(DNA_Utility::map_base(s[kmer_idx + k]));
        kmer_idx++;
    }
}


template <uint16_t k, uint16_t l>
inline bool Streaming_Kmer_Hashset<k, l>::find(const Kmer<k> key, const std::size_t b) const
{
    assert(b < cap_);   assert(S[b] <= B);

    for(std::size_t j = 0; j < S[b]; ++j)
        if(K[b * B + j] == key)
            return true;

    return false;
}


template <uint16_t k, uint16_t l>
inline bool Streaming_Kmer_Hashset<k, l>::insert_at(const Kmer<k> key, const std::size_t b)
{
    if(find(key, b))
        return true;

    if(S[b] == B)
        return false;

    K[b * B + S[b]] = key;
    S[b]++;
    return true;
}


template <uint16_t k, uint16_t l>
inline void Streaming_Kmer_Hashset<k, l>::preprocess_skew_mins(const std::vector<std::string>& file_paths)
{
    std::size_t scanned_kmer = 0;
    constexpr std::size_t scan_th = 100'000;

    Min_Iterator<k, l> min_it;
    ankerl::unordered_dense::map<uint64_t, uint64_t> min_freq;

    for(const auto& path : file_paths)
    {
        Parser input(path);
        while(scanned_kmer < scan_th && input.read_next_seq())
        {
            auto const seq = input.seq();
            auto const len = input.seq_len();
            if(len < k)
                continue;

            min_it.reset(seq);
            std::size_t kmer_idx = 0;
            while(scanned_kmer < scan_th)
            {
                min_freq[min_it.hash()]++;

                if(kmer_idx + k == len)
                    break;

                min_it.advance(seq[kmer_idx + k]);
                kmer_idx++;
                scanned_kmer++;
            }
        }
    }


    for(auto [min, f] : min_freq)
        if(f > skew_th)
            skew_mins.insert(min);

    std::cerr << "#skew-minimizers: " << skew_mins.size() << "\n";
}


template <uint16_t k, uint16_t l>
inline void Streaming_Kmer_Hashset<k, l>::report_stats() const
{
    std::cerr   << "Table height: " << cap_ << ".\n";
    std::cerr   << "Table width:  " << B << ".\n";
    std::cerr   << "========================================\n";

    std::cerr   << "Total sequence length:  " << sum_len << ".\n";
    std::cerr   << "Number of super k-mers: " << sup_kmer_c << ".\n";
    std::cerr   << "========================================\n";


    uint64_t sz = 0;
    std::for_each(S, S + cap_, [&](auto s){ sz += s; });

    uint64_t skew_min = 0;
    for(const auto& m : min_to_kmer)
        if(m.second.size() > B)
            skew_min++;

    std::cerr   << "Minimizer set size: " << min_to_kmer.size() << ".\n";
    std::cerr   << "Minimizers with more than " << B << " k-mers: " << skew_min << ".\n";
    std::cerr   << "========================================\n";

    std::cerr   << "Complete set size: " << sz << ".\n";
    std::cerr   << "Overflow set size: " << overflow.size() << " ("
                    << (overflow.size() * 100.0 / (sz + overflow.size())) << "%).\n";
    std::cerr   << "========================================\n";


    uint64_t full = 0, empty = 0;
    for(std::size_t b = 0; b < cap_; ++b)
        full += (S[b] == B), empty += (S[b] == 0);

    std::cerr   << "#Blocks:       " << cap_ << ".\n";
    std::cerr   << "#Full-blocks:  " << full << "(" << (full * 100.0 / cap_) << "%).\n";
    std::cerr   << "#Empty-blocks: " << empty << "(" << (empty * 100.0 / cap_) << "%).\n";
    std::cerr   << "========================================\n";

}



#endif
