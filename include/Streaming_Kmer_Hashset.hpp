
#ifndef STREAMING_KMER_HASHSET_HPP
#define STREAMING_KMER_HASHSET_HPP



#include "DNA_Utility.hpp"
#include "Kmer.hpp"
#include "Directed_Vertex.hpp"
#include "Rolling_Hash.hpp"
#include "Minimizer_Iterator.hpp"
#include "Hash_Table.hpp"
#include "utility.hpp"

#include "unordered_dense/unordered_dense.h"

#include <immintrin.h>

#include <cstdint>
#include <cstddef>
#include <bit>
#include <cstring>
#include <iostream>
#include <cassert>


// Hashset particularly suited for streaming k-mer operations.
template <uint16_t k, uint16_t l = 19>
class Streaming_Kmer_Hashset
{
private:

    static constexpr double lf_default = 0.8;   // Default load factor: 80%.
    static constexpr std::size_t B = 32;    // Bucket size.
    static_assert(B % 8 == 0, "Bucket size must be a multiple of 8.");
    static_assert(B == 32, "Bucket size needs to be 32 for the current vectorization scheme to function.");

    std::size_t cap_;   // Capacity of the table in terms of buckets; adjusted to be a power of 2.

    Kmer<k>* K; // Flat table of key-collection of size `cap x B`.

    typedef uint8_t checksum_t;
    checksum_t* C;  // Checksum table for the keys.

    uint64_t sz;    // Size of the main table.

    static constexpr uint64_t skew_th = B;
    Hash_Set<Kmer<k>> overflow; // The overflow table.

    uint64_t sum_len;   // Total sequence length provided.
    uint64_t sup_kmer_c;    // Number of super k-mers in the provided sequences to the set.

    // Following are for used for developer stats.
    ankerl::unordered_dense::map<uint64_t, ankerl::unordered_dense::set<Kmer<k>, Kmer_Hasher<k>>> min_to_kmer;  // Set of overflowing minimizers.


    // Returns the size of the `b`'th bucket.
    std::size_t bucket_size(std::size_t b) const;

    // Returns `true` iff the k-mer `key` with checksum `c` is present at the
    // `i`'th bucket.
    bool find(Kmer<k> key, checksum_t c, std::size_t i) const;

    // Directly inserts the k-mer `key` with checksum `c` into the `i`'th
    // bucket, at index `j` of the bucket.
    void insert_at(const Kmer<k> key, const checksum_t c, const std::size_t i, const std::size_t j) { K[i * B + j] = key, C[i * B + j] = c, sz++; }

    // Tries to insert the k-mer `key` with checksum `c` into the `i`'th bucket,
    // returning `true` iff `key` is absent at the bucket after this operation.
    bool try_insert_at(Kmer<k> key, checksum_t c, std::size_t i);

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
    auto size() const { return sz + overflow_size(); }

    // Returns the size of the overflow table.
    auto overflow_size() const { return overflow.size(); }

    // Clears the table.
    void clear();

    // Tries to insert all the k-mers from the `l`-length sequence `S` into the
    // table.
    template <bool track_stats_ = false> void insert(const char* s, std::size_t len);

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
    , C(allocate<checksum_t>(cap_ * B))
    , sz(0)
{
    clear();
}


template <uint16_t k, uint16_t l>
inline Streaming_Kmer_Hashset<k, l>::~Streaming_Kmer_Hashset()
{
    deallocate(K), deallocate(C);
}


template <uint16_t k, uint16_t l>
inline void Streaming_Kmer_Hashset<k, l>::clear()
{
    sum_len = 0, sup_kmer_c = 0;

    std::memset(C, 0, cap_ * B);    // Checksum `0` denotes absence of a key.
    sz = 0;

    // skew_mins.clear();
    overflow.clear();

    min_to_kmer.clear();
}


template <uint16_t k, uint16_t l>
inline bool Streaming_Kmer_Hashset<k, l>::find(const Kmer<k> key, const uint8_t c, const std::size_t i) const
{
    assert(i < cap_);

    const auto c_vec = _mm256_load_si256(reinterpret_cast<const __m256i*>(C + (i * B)));
    const auto query = _mm256_set1_epi8(static_cast<char>(c));
    auto cmp = _mm256_cmpeq_epi8_mask(c_vec, query);

    while(cmp)
    {
        const auto j = __builtin_ctz(cmp);
        if(K[i * B + j] == key)
            return true;

        cmp &= (cmp - 1);
    }

    return false;
}


template <uint16_t k, uint16_t l>
inline bool Streaming_Kmer_Hashset<k, l>::try_insert_at(const Kmer<k> key, const uint8_t c, const std::size_t i)
{
    if(find(key, c, i))
        return true;

    const auto c_vec = _mm256_load_si256(reinterpret_cast<const __m256i*>(C + (i * B)));
    const auto zero = _mm256_setzero_si256();
    const auto cmp = _mm256_cmpeq_epi8_mask(c_vec, zero);

    const auto j = std::countr_zero(cmp);
    if(j == 32)
        return false;

    insert_at(key, c, i, j);
    return true;
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
    Rolling_Hash<k, true> rh(s);
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
        const auto c = std::max(rh.template checksum<sizeof(checksum_t) * 8>(), uint64_t(1));   // Avoiding checksum 0 by overloading checksum 1.
        if(!try_insert_at(key, c, b))
        {
            // TODO: get the following hashes from the minimizer iterator (rolling hashes), instead of full-blown xxHash computations.
            const auto h_1 = XXH3_64bits_withSeed(&h, sizeof(h), 0);
            const auto h_2 = XXH3_64bits_withSeed(&h, sizeof(h), 1llu << 63);

            const auto p = h_1 & idx_mask;
            const auto q = h_2 & idx_mask;
            if(!find(key, c, p) && !find(key, c, q))
            {
                const auto zero = _mm256_setzero_si256();

                const auto c_vec_p = _mm256_load_si256(reinterpret_cast<const __m256i*>(C + (p * B)));
                const auto c_vec_q = _mm256_load_si256(reinterpret_cast<const __m256i*>(C + (q * B)));

                const auto cmp_p = _mm256_cmpeq_epi8_mask(c_vec_p, zero);
                const auto cmp_q = _mm256_cmpeq_epi8_mask(c_vec_q, zero);

                const auto sz_p = B - std::popcount(cmp_p);
                const auto sz_q = B - std::popcount(cmp_q);

                if(sz_p < sz_q)
                    insert_at(key, c, p, std::countr_zero(cmp_p));
                else if(sz_q < B)
                    insert_at(key, c, q, std::countr_zero(cmp_q));
                else
                    overflow.insert(key);
            }
        }

        if(kmer_idx + k == len)
            break;

        v.roll_forward(DNA_Utility::map_base(s[kmer_idx + k]));
        min_it.advance(s[kmer_idx + k]);
        rh.advance(s[kmer_idx + k]);
        kmer_idx++;
    }
}


template <uint16_t k, uint16_t l>
inline std::size_t Streaming_Kmer_Hashset<k, l>::bucket_size(const std::size_t b) const
{
    const auto zero  = _mm256_setzero_si256();
    const auto c_vec = _mm256_load_si256(reinterpret_cast<const __m256i*>(C + (b * B)));
    const auto cmp   = _mm256_cmpeq_epi8_mask(c_vec, zero);
    return B - std::popcount(cmp);
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


    uint64_t sz = size();

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
        full += (bucket_size(b) == B), empty += (bucket_size(b) == 0);

    std::cerr   << "#Blocks:       " << cap_ << ".\n";
    std::cerr   << "#Full-blocks:  " << full << "(" << (full * 100.0 / cap_) << "%).\n";
    std::cerr   << "#Empty-blocks: " << empty << "(" << (empty * 100.0 / cap_) << "%).\n";
    std::cerr   << "========================================\n";

}



#endif
