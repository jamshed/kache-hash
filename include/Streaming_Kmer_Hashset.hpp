
#ifndef STREAMING_KMER_HASHSET_HPP
#define STREAMING_KMER_HASHSET_HPP



#include "DNA_Utility.hpp"
#include "Kmer.hpp"
#include "Directed_Vertex.hpp"
#include "Minimizer_Iterator.hpp"
#include "Hash_Table.hpp"
#include "utility.hpp"

#include "unordered_dense/unordered_dense.h"
#include "boost/unordered/concurrent_flat_set.hpp"

#include <immintrin.h>

#include <cstdint>
#include <cstddef>
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
    static_assert(B % 8 == 0, "Bucket size must be a multiple of 8.");

    std::size_t cap_;   // Capacity of the table in terms of buckets; adjusted to be a power of 2.

    Kmer<k>* K; // Flat table of key-collection of size `cap x B`.
    uint8_t* S; // Size of each bucket.

    static constexpr uint64_t skew_th = B;
    // Hash_Set<Kmer<k>> overflow; // The overflow table.
    boost::unordered::concurrent_flat_set<Kmer<k>, Hash<Kmer<k>>> overflow;

    uint64_t sum_len;   // Total sequence length provided.
    uint64_t sup_kmer_c;    // Number of super k-mers in the provided sequences to the set.

    static constexpr uint8_t lock_mask = 0b1000'0000;   // Bit of size-counter used to lock the bucket.

    // Following are for used for developer stats.
    ankerl::unordered_dense::map<uint64_t, ankerl::unordered_dense::set<Kmer<k>, Kmer_Hasher<k>>> min_to_kmer;  // Set of overflowing minimizers.


    // Returns the size of the `b`'th bucket.
    std::size_t size(const std::size_t b) const { return S[b] & ~lock_mask; }

    // Tries to insert the k-mer `key` into the bucket `b`, returning `true` iff
    // `key` is absent at bucket `b` after this operation.
    bool insert_at(Kmer<k> key, std::size_t b);

    // Returns `true` iff the k-mer `key` is present at the bucket `b`.
    bool find(Kmer<k> key, std::size_t b) const;

    // Returns `true` iff the k-mer `key` is present at the bucket `b`.
    bool find_512(Kmer<k> key, std::size_t b) const;

    // Locks the `b`'th bucket for exclusive write-access.
    void lock(std::size_t b);

    // Unlocks the `b`'th bucket from exclusive write-access.
    void unlock(std::size_t b);

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

    // Returns the size of the overflow table.
    auto overflow_size() const { return overflow.size(); }

    // Clears the table.
    void clear();

    // Tries to insert all the k-mers from the `l`-length sequence `S` into the
    // table.
    template <bool track_stats_ = false> void insert(const char* s, std::size_t len);

    // Returns the number of k-mer instances from the `l`-length sequence `s`
    // found in the table.
    std::size_t find(const char* s, std::size_t len) const;

/*
    // Preprocesses a prefix of sequences from the files at `file_path` to
    // populate a collection of skew-minimizers.
    void preprocess_skew_mins(const std::vector<std::string>& file_paths);

    // Preprocesses a random set of k-mers from the files at `file_path` to
    // populate a collection of skew-minimizers. Proof-of-concept method.
    void preprocess_skew_mins_random(const std::vector<std::string>& file_paths);
*/

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

    // skew_mins.clear();
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

        const auto b = h & idx_mask;
        if(!insert_at(key, b))
        {
            const auto h_1 = XXH3_64bits_withSeed(&h, sizeof(h), 0);
            const auto h_2 = XXH3_64bits_withSeed(&h, sizeof(h), 1llu << 63);

            const auto s_1 = h_1 & idx_mask;
            const auto s_2 = h_2 & idx_mask;
            assert(s_1 != s_2);
            const auto p = std::min(s_1, s_2), q = std::max(s_1, s_2);
            lock(p), lock(q);
            const auto s_p = size(p), s_q = size(q);
            if(find(key, p) || find(key, q))
                unlock(q), unlock(p);
            else
            {
                if(s_p < s_q)
                {
                    K[p * B + s_p] = key, S[p]++;
                    unlock(q), unlock(p);
                }
                else if(s_q < B)
                {
                    K[q * B + s_q] = key, S[q]++;
                    unlock(q), unlock(p);
                }
                else
                {
                    unlock(q), unlock(p);
                    overflow.insert(key);   // TODO: fix.
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
inline std::size_t Streaming_Kmer_Hashset<k, l>::find(const char* const s, const std::size_t len) const
{
    Directed_Vertex<k> v(s);
    std::size_t kmer_idx = 0;
    std::size_t c = 0;

    Min_Iterator<k, l> min_it;
    min_it.reset(s);
    while(true)
    {
        const auto key = v.canonical();
        const auto h = min_it.hash();
        const auto idx_mask = cap_ - 1;

        const auto b = h & idx_mask;
        if(find(key, b))
            c++;
        else if(size(b) == B)
        {
            const auto p = XXH3_64bits_withSeed(&h, sizeof(h), 0) & idx_mask;
            const auto q = XXH3_64bits_withSeed(&h, sizeof(h), 1llu << 63) & idx_mask;
            if(find(key, p) || find(key, q))
                c++;
            else if(size(p) == B && size(q) == B)
                c += overflow.count(key);
        }

        if(kmer_idx + k == len)
            break;

        min_it.advance(s[kmer_idx + k]);
        v.roll_forward(DNA_Utility::map_base(s[kmer_idx + k]));
        kmer_idx++;
    }

    return c;
}


template <uint16_t k, uint16_t l>
inline bool Streaming_Kmer_Hashset<k, l>::find(const Kmer<k> key, const std::size_t b) const
{
/*
    assert(b < cap_);   assert(S[b] <= B);

    for(std::size_t j = 0; j < S[b]; ++j)
        if(K[b * B + j] == key)
            return true;

    return false;
*/

    return find_512(key, b);
}


template <uint16_t k, uint16_t l>
inline bool Streaming_Kmer_Hashset<k, l>::find_512(const Kmer<k> key, const std::size_t b) const
{
    assert(b < cap_);   assert(size(b) <= B);

    const auto query = _mm512_set1_epi64(key.as_int());

    std::size_t j = 0;
    const auto sz = size(b);
    for(; j < sz && j + 8 <= B; j += 8)
    {
        const auto v = _mm512_load_si512(static_cast<const void*>(K + (b * B + j)));
        const auto cmp = _mm512_cmpeq_epi64_mask(v, query);
        if(cmp)
        {
            const auto off = __builtin_ctz(static_cast<uint8_t>(cmp));
            return static_cast<std::size_t>(off) < sz;
        }
    }

    return false;
}


template <uint16_t k, uint16_t l>
inline void Streaming_Kmer_Hashset<k, l>::lock(const std::size_t b)
{
    while(true)
    {
        const uint8_t exp_unlocked = S[b] & ~lock_mask;
        if(__sync_bool_compare_and_swap(S + b, exp_unlocked, exp_unlocked | lock_mask))
            break;
    }
}


template <uint16_t k, uint16_t l>
inline void Streaming_Kmer_Hashset<k, l>::unlock(const std::size_t b)
{
    assert(S[b] & lock_mask);
    S[b] &= ~lock_mask;
}


template <uint16_t k, uint16_t l>
inline bool Streaming_Kmer_Hashset<k, l>::insert_at(const Kmer<k> key, const std::size_t b)
{
    lock(b);

    if(find(key, b))
    {
        unlock(b);
        return true;
    }

    const auto sz = size(b);
    if(sz == B)
    {
        unlock(b);
        return false;
    }

    K[b * B + sz] = key;
    S[b]++;
    unlock(b);
    return true;
}


/*
template <uint16_t k, uint16_t l>
inline void Streaming_Kmer_Hashset<k, l>::preprocess_skew_mins(const std::vector<std::string>& file_paths)
{
    std::size_t scanned_kmer = 0;
    constexpr std::size_t scan_th = 1'000'000;

    Min_Iterator<k, l> min_it;
    Hash_Map<uint64_t, uint64_t> min_freq;

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
    std::cerr << "Minimizer-counting map size: " << min_freq.size() << "\n";
    std::cerr << "Minimizer-counting map cap:  " << min_freq.capacity() << "\n";
}


template <uint16_t k, uint16_t l>
inline void Streaming_Kmer_Hashset<k, l>::preprocess_skew_mins_random(const std::vector<std::string>& file_paths)
{
    constexpr std::size_t scan_th = 1'000'000;

    std::string s;
    for(const auto& path : file_paths)
    {
        Parser input(path);
        while(input.read_next_seq())
            s += input.seq();
    }

    Min_Iterator<k, l> min_it;
    Hash_Map<uint64_t, uint64_t> min_freq;

    std::mt19937 r(std::chrono::system_clock::now().time_since_epoch().count());
    for(std::size_t i = 0; i < scan_th; ++i)
    {
        const auto idx = r() % (s.length() - k);
        min_it.reset(s.c_str() + idx);
        min_freq[min_it.hash()]++;
    }


    for(auto [min, f] : min_freq)
        if(f > skew_th)
            skew_mins.insert(min);

    std::cerr << "#skew-minimizers: " << skew_mins.size() << "\n";
    std::cerr << "Minimizer-counting map size: " << min_freq.size() << "\n";
    std::cerr << "Minimizer-counting map cap:  " << min_freq.capacity() << "\n";
}
*/


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
    uint64_t skew_kmer = 0;
    for(const auto& m : min_to_kmer)
        if(m.second.size() > skew_th)
            skew_min++, skew_kmer += m.second.size();

    std::cerr   << "Minimizer set size: " << min_to_kmer.size() << ".\n";
    std::cerr   << "Minimizers with more than " << B << " k-mers: " << skew_min << " (" << (skew_min * 100.0 / min_to_kmer.size()) << "%).\n";
    std::cerr   << "Associated k-mer count: " << skew_kmer << " (" << (skew_kmer * 100.0 / sz) << "%).\n";
    // std::cerr   << "Identified skew-minimizer set size: " << skew_mins.size() << "\n";
    // std::cerr   << "Identified skew-minimizer data structure size: " << (skew_mins.capacity() * sizeof(uint64_t) / (1024lu))  << " KBs.\n";
    std::cerr   << "========================================\n";

    std::cerr   << "Complete set size: " << sz << ".\n";
    std::cerr   << "Overflow set size: " << overflow.size() << " ("
                    << (overflow.size() * 100.0 / (sz + overflow.size())) << "%).\n";
    std::cerr   << "========================================\n";


    uint64_t full = 0, empty = 0;
    for(std::size_t b = 0; b < cap_; ++b)
        full += (size(b) == B), empty += (size(b) == 0);

    std::cerr   << "#Blocks:       " << cap_ << ".\n";
    std::cerr   << "#Full-blocks:  " << full << "(" << (full * 100.0 / cap_) << "%).\n";
    std::cerr   << "#Empty-blocks: " << empty << "(" << (empty * 100.0 / cap_) << "%).\n";
    std::cerr   << "========================================\n";

}



#endif
