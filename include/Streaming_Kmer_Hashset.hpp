
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
#include "boost/unordered/concurrent_flat_set.hpp"

#include <immintrin.h>

#include <cstdint>
#include <cstddef>
#include <bit>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>


// Hashset particularly suited for streaming k-mer operations.
// NB: concurrent insertions are guaranteed to function correctly only if
// deletions are not performed.
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

    // Metadata of the k-mers in a bucket.
    struct Metadata
    {
        uint8_t cs[B];  // Checksums of the k-mers.
        uint8_t min_coord[B];   // Minimizer coordinates within the k-mers.

        static constexpr uint8_t lock_mask = 0b1000'0000;   // Bit of `min_coord[0]` used to lock the bucket.
    };

    Metadata* M;    // Metadata table of the buckets.

    static constexpr uint64_t skew_th = B;
    Hash_Set<Kmer<k>> overflow; // The overflow table.
    boost::unordered::concurrent_flat_set<Kmer<k>, Hash<Kmer<k>>> overflow_mt;  // The overflow table in multithreaded setting.

    uint64_t sum_len;   // Total sequence length provided.
    uint64_t sup_kmer_c;    // Number of super k-mers in the provided sequences to the set.

    // Following are for used for developer stats.
    ankerl::unordered_dense::map<uint64_t, ankerl::unordered_dense::set<Kmer<k>, Kmer_Hasher<k>>> min_to_kmer;  // Set of overflowing minimizers.


    // Returns the size of the `b`'th bucket.
    std::size_t bucket_size(std::size_t b) const;

    // Locks the `b`'th bucket for exclusive-write access.
    void lock(std::size_t b);

    // Unlocks the `b`'th bucket from exclusive-write access.
    void unlock(std::size_t b);

    // Returns `true` iff the k-mer `key` with checksum `c` is present in the
    // `i`'th bucket, whose checksum-vector is `c_vec`.
    bool find(Kmer<k> key, uint8_t c, std::size_t i, __m256i c_vec) const;

    // Returns `true` iff the k-mer `key` with checksum `c` is present in the
    // `i`'th bucket.
    bool find(Kmer<k> key, uint8_t c, std::size_t i) const;

    // Directly places the k-mer `key` with checksum `c` into the `i`'th
    // bucket at index `j`.
    void place_at(const Kmer<k> key, const uint8_t c, const std::size_t i, const std::size_t j) { K[i * B + j] = key, M[i].cs[j] = c; }

    // Tries to insert the k-mer `key` with checksum `c` into the `i`'th bucket,
    // returning `true` iff `key` is present at the bucket after this operation.
    // `mt_` denotes whether multiple threads can access the table concurrently.
    template <bool mt_> bool try_insert_at(Kmer<k> key, uint8_t c, std::size_t i);

public:

    // Constructs a streaming k-mer hashset to support up-to `max_sz` k-mers,
    // with a maximum load-factor of `lf`.
    Streaming_Kmer_Hashset(std::size_t max_sz, double lf = lf_default);

    // Deserializes a streaming k-mer hashset from the stream `is`.
    Streaming_Kmer_Hashset(std::ifstream& is);

    Streaming_Kmer_Hashset(const Streaming_Kmer_Hashset&) = delete;
    Streaming_Kmer_Hashset(Streaming_Kmer_Hashset&&) = delete;
    Streaming_Kmer_Hashset& operator=(const Streaming_Kmer_Hashset&) = delete;
    Streaming_Kmer_Hashset& operator=(Streaming_Kmer_Hashset&&) = delete;

    ~Streaming_Kmer_Hashset();

    // Returns the capacity (in terms of keys) of the table.
    auto capacity() const { return cap_ * B; }

    // Returns the size of the table. `mt_` denotes whether multiple threads can
    // access the table concurrently.
    template <bool mt_> std::size_t size() const;

    // Returns the size of the overflow table.
    template <bool mt_> auto overflow_size() const { return !mt_ ? overflow.size() : overflow_mt.size();  }

    // Clears the table.
    void clear();

    // Tries to insert all the k-mers from the `len`-length sequence `s` into
    // the table.
    template <bool track_stats_ = false> void insert(const char* s, std::size_t len);

    // Tries to insert all the k-mers from the `len`-length sequence `s` into
    // the table in a concurrent manner.
    void insert_mt(const char* s, std::size_t len);

    // Returns the number of k-mer instances from the `len`-length sequence `s`
    // present in the set.
    std::size_t find(const char* s, std::size_t len) const;

    // Serializes the table to the stream `os`.
    void serialize(std::ofstream& os) const;

    // Deserializes the table from the stream `is`.
    void deserialize(std::ifstream& is);

    // Reports developer statistics about the table.
    void report_stats() const;
};


template <uint16_t k, uint16_t l>
inline Streaming_Kmer_Hashset<k, l>::Streaming_Kmer_Hashset(const std::size_t max_sz, const double lf):
      cap_(ceil_pow_2((max_sz / lf) / B))
    , K(aligned_allocate<Kmer<k>>(cap_ * B, L1_CACHE_LINE_SIZE))
    , M(aligned_allocate<Metadata>(cap_, L1_CACHE_LINE_SIZE))
{
    clear();
}


template <uint16_t k, uint16_t l>
inline Streaming_Kmer_Hashset<k, l>::Streaming_Kmer_Hashset(std::ifstream& is):
      cap_(0)
    , K(nullptr)
    , M(nullptr)
{
    deserialize(is);
}


template <uint16_t k, uint16_t l>
inline Streaming_Kmer_Hashset<k, l>::~Streaming_Kmer_Hashset()
{
    deallocate(K), deallocate(M);
}


template <uint16_t k, uint16_t l>
inline void Streaming_Kmer_Hashset<k, l>::clear()
{
    sum_len = 0, sup_kmer_c = 0;

    std::memset(M, 0, cap_ * sizeof(Metadata)); // Checksum `0` denotes absence of a key.

    overflow.clear(), overflow_mt.clear();

    min_to_kmer.clear();
}


template <uint16_t k, uint16_t l>
inline std::size_t Streaming_Kmer_Hashset<k, l>::bucket_size(const std::size_t b) const
{
    const auto zero  = _mm256_setzero_si256();
    const auto c_vec = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[b].cs));
    const auto cmp   = _mm256_cmpeq_epi8_mask(c_vec, zero);

    assert(std::countr_zero(cmp) + std::countl_one(cmp) == B);
    return std::countr_zero(cmp);   // This must hold if no deletion is performed at the table.
}


template <uint16_t k, uint16_t l>
inline void Streaming_Kmer_Hashset<k, l>::lock(const std::size_t b)
{
    while(true)
    {
        const auto exp_unlocked = M[b].min_coord[0] & ~Metadata::lock_mask;
        if(__sync_bool_compare_and_swap(&M[b].min_coord[0], exp_unlocked, exp_unlocked | Metadata::lock_mask))
            break;
    }
}


template <uint16_t k, uint16_t l>
inline void Streaming_Kmer_Hashset<k, l>::unlock(const std::size_t b)
{
    assert(M[b].min_coord[0] & Metadata::lock_mask);
    M[b].min_coord[0] &= ~Metadata::lock_mask;
}


template <uint16_t k, uint16_t l>
template <bool mt_>
inline std::size_t Streaming_Kmer_Hashset<k, l>::size() const
{
    std::size_t sz = 0;
    for(std::size_t b = 0; b < cap_; ++b)
        sz += bucket_size(b);

    return sz + overflow_size<mt_>();
}


template <uint16_t k, uint16_t l>
inline bool Streaming_Kmer_Hashset<k, l>::find(const Kmer<k> key, const uint8_t c, const std::size_t i, const __m256i c_vec) const
{
    assert(i < cap_);

    const auto query = _mm256_set1_epi8(static_cast<char>(c));
    auto cmp = _mm256_cmpeq_epi8_mask(c_vec, query);

    while(cmp)
    {
        const auto j = std::countr_zero(cmp);
        if(K[i * B + j] == key)
            return true;

        cmp &= (cmp - 1);
    }

    return false;
}


template <uint16_t k, uint16_t l>
inline bool Streaming_Kmer_Hashset<k, l>::find(const Kmer<k> key, const uint8_t c, const std::size_t i) const
{
    const auto c_vec = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[i].cs));
    return find(key, c, i, c_vec);
}


template <uint16_t k, uint16_t l>
template <bool mt_>
inline bool Streaming_Kmer_Hashset<k, l>::try_insert_at(const Kmer<k> key, const uint8_t c, const std::size_t i)
{
    if constexpr(mt_)   lock(i);
    const auto c_vec = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[i].cs));

    if(find(key, c, i, c_vec))
    {
        if constexpr(mt_)   unlock(i);
        return true;
    }

    const auto zero = _mm256_setzero_si256();
    const auto cmp = _mm256_cmpeq_epi8_mask(c_vec, zero);

    const auto j = std::countr_zero(cmp);
    if(j == B)
    {
        if constexpr(mt_)   unlock(i);
        return false;
    }

    place_at(key, c, i, j);
    if constexpr(mt_)   unlock(i);
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
        const auto c = std::max(rh.template checksum<8>(), uint64_t(1));    // Avoiding checksum 0 by overloading checksum 1.
        if(!try_insert_at<false>(key, c, b))
        {
            // TODO: get the following hashes from the minimizer iterator (rolling hashes), instead of full-blown xxHash computations.
            const auto h_1 = XXH3_64bits_withSeed(&h, sizeof(h), 0);
            const auto h_2 = XXH3_64bits_withSeed(&h, sizeof(h), 1llu << 63);

            const auto p = h_1 & idx_mask;
            const auto q = h_2 & idx_mask;
            if(!find(key, c, p) && !find(key, c, q))
            {
                const auto zero = _mm256_setzero_si256();

                const auto c_vec_p = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[p].cs));
                const auto c_vec_q = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[q].cs));

                const auto cmp_p = _mm256_cmpeq_epi8_mask(c_vec_p, zero);
                const auto cmp_q = _mm256_cmpeq_epi8_mask(c_vec_q, zero);

                const auto sz_p = static_cast<std::size_t>(std::countr_zero(cmp_p));
                const auto sz_q = static_cast<std::size_t>(std::countr_zero(cmp_q));
                assert(std::countr_zero(cmp_p) + std::countl_one(cmp_p) == B), assert(std::countr_zero(cmp_q) + std::countl_one(cmp_q) == B);

                if(sz_p < sz_q)
                    place_at(key, c, p, sz_p);
                else if(sz_q < B)
                    place_at(key, c, q, sz_q);
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
inline void Streaming_Kmer_Hashset<k, l>::insert_mt(const char* const s, const std::size_t len)
{
    // sum_len += len;
    if(len < k)
        return;

    Directed_Vertex<k> v(s);
    std::size_t kmer_idx = 0;

    Min_Iterator<k, l> min_it;
    Rolling_Hash<k, true> rh(s);
    min_it.reset(s);
    while(true)
    {
        const auto key = v.canonical();
        const auto h = min_it.hash();
        const auto idx_mask = cap_ - 1;

        const auto b = h & idx_mask;
        const auto c = std::max(rh.template checksum<8>(), uint64_t(1));    // Avoiding checksum 0 by overloading checksum 1.
        if(!try_insert_at<true>(key, c, b))
        {
            // TODO: get the following hashes from the minimizer iterator (rolling hashes), instead of full-blown xxHash computations.
            const auto h_1 = XXH3_64bits_withSeed(&h, sizeof(h), 0);
            const auto h_2 = XXH3_64bits_withSeed(&h, sizeof(h), 1llu << 63);

            const auto s_1 = h_1 & idx_mask;
            const auto s_2 = h_2 & idx_mask;
            const auto p = std::min(s_1, s_2), q = std::max(s_1, s_2);
            assert(p != q);

            lock(p), lock(q);
            while(true)
                if(find(key, c, p) || find(key, c, q))
                {
                    unlock(q), unlock(p);
                    break;
                }
                else
                {
                    const auto zero = _mm256_setzero_si256();

                    const auto c_vec_p = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[p].cs));
                    const auto c_vec_q = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[q].cs));

                    const auto cmp_p = _mm256_cmpeq_epi8_mask(c_vec_p, zero);
                    const auto cmp_q = _mm256_cmpeq_epi8_mask(c_vec_q, zero);

                    const auto sz_p = static_cast<std::size_t>(std::countr_zero(cmp_p));
                    const auto sz_q = static_cast<std::size_t>(std::countr_zero(cmp_q));
                    assert(std::countr_zero(cmp_p) + std::countl_one(cmp_p) == B), assert(std::countr_zero(cmp_q) + std::countl_one(cmp_q) == B);

                    if(sz_p < sz_q)
                    {
                        place_at(key, c, p, sz_p);
                        unlock(q), unlock(p);
                        break;
                    }
                    else if(sz_q < B)
                    {
                        place_at(key, c, q, sz_q);
                        unlock(q), unlock(p);
                        break;
                    }
                    else
                    {
                        unlock(q), unlock(p);
                        overflow_mt.insert(key);
                        break;
                    }
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
inline std::size_t  Streaming_Kmer_Hashset<k, l>::find(const char* const s, const std::size_t len) const
{
    Directed_Vertex<k> v(s);
    std::size_t kmer_idx = 0;
    std::size_t pos_q = 0;

    Min_Iterator<k, l> min_it;
    Rolling_Hash<k, true> rh(s);
    min_it.reset(s);
    while(true)
    {
        const auto key = v.canonical();
        const auto h = min_it.hash();
        const auto idx_mask = cap_ - 1;

        const auto b = h & idx_mask;
        const auto c = std::max(rh.template checksum<8>(), uint64_t(1));    // Avoiding checksum 0 by overloading checksum 1.
        if(find(key, c, b))
            pos_q++;
        else
        {
            // TODO: get the following hashes from the minimizer iterator (rolling hashes), instead of full-blown xxHash computations.
            const auto h_1 = XXH3_64bits_withSeed(&h, sizeof(h), 0);
            const auto h_2 = XXH3_64bits_withSeed(&h, sizeof(h), 1llu << 63);

            const auto p = h_1 & idx_mask;
            const auto q = h_2 & idx_mask;
            if(find(key, c, p) || find(key, c, q))
                pos_q++;
            else if(bucket_size(p) == B && bucket_size(q) == B)
                pos_q += overflow.contains(key);
        }

        if(kmer_idx + k == len)
            break;

        v.roll_forward(DNA_Utility::map_base(s[kmer_idx + k]));
        min_it.advance(s[kmer_idx + k]);
        rh.advance(s[kmer_idx + k]);
        kmer_idx++;
    }

    return pos_q;
}


template <uint16_t k, uint16_t l>
inline void Streaming_Kmer_Hashset<k, l>::serialize(std::ofstream& os) const
{
    ::serialize(k, os), ::serialize(l, os), ::serialize(B, os);

    ::serialize(cap_, os);
    ::serialize(K, cap_ * B, os), ::serialize(M, cap_, os);

    overflow.serialize(os);

    ::serialize(sum_len, os), ::serialize(sup_kmer_c, os);
}


template <uint16_t k, uint16_t l>
inline void Streaming_Kmer_Hashset<k, l>::deserialize(std::ifstream& is)
{
    uint16_t k_, l_;
    std::size_t B_;
    ::deserialize(k_, is), ::deserialize(l_, is), ::deserialize(B_, is);
    if(k_ != k || l_ != l || B_ != B)
    {
        std::cerr << "Mismatching configuration of the hash table. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    ::deserialize(cap_, is);

    if(K)
        deallocate(K), deallocate(M);

    K = aligned_allocate<Kmer<k>>(cap_ * B, L1_CACHE_LINE_SIZE), M = aligned_allocate<Metadata>(cap_, L1_CACHE_LINE_SIZE);
    ::deserialize(K, cap_ * B, is), ::deserialize(M, cap_, is);

    overflow.deserialize(is);

    ::deserialize(sum_len, is), ::deserialize(sup_kmer_c, is);
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


    uint64_t sz = size<false>();

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
                    << (overflow.size() * 100.0 / sz) << "%).\n";
    std::cerr   << "========================================\n";


    uint64_t full = 0, empty = 0, blk_occ = 0;
    for(std::size_t b = 0; b < cap_; ++b)
        full += (bucket_size(b) == B), empty += (bucket_size(b) == 0), blk_occ += bucket_size(b);

    std::cerr   << "#Blocks:       " << cap_ << ".\n";
    std::cerr   << "#Full-blocks:  " << full << "(" << (full * 100.0 / cap_) << "%).\n";
    std::cerr   << "#Empty-blocks: " << empty << "(" << (empty * 100.0 / cap_) << "%).\n";
    std::cerr   << "Avg-occupancy: " << blk_occ / cap_ << ".\n";
    std::cerr   << "========================================\n";

}



#endif
