
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
#include "boost/unordered/concurrent_flat_map.hpp"

#include <cstdint>
#include <cstddef>
#include <bit>
#include <optional>
#include <utility>
#include <type_traits>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>

#include <immintrin.h>


template <uint16_t k, uint16_t l> class Kmer_Window;


// Hash table particularly suited for streaming `k`-mer operations. `mt_`
// denotes whether a concurrent table is desired. If a hash map is required, the
// mapped type is `T_`. `l`-minimizers are used in underlying computations.
// NB: concurrent insertions are guaranteed to function correctly only if
// deletions are not performed, as of now.
template <uint16_t k, bool mt_, typename T_ = void, uint16_t l = 19>
class Streaming_Kmer_Hash_Table
{
private:

    static constexpr double lf_default = 0.8;   // Default load factor: 80%.
    static constexpr std::size_t B = 32;    // Bucket size.
    static_assert(B % 8 == 0, "Bucket size must be a multiple of 8.");
    static_assert(B == 32, "Bucket size needs to be 32 for the current vectorization scheme to function.");

    static constexpr bool is_map_ = !std::is_void_v<T_>;    // Whether a hash-map or a -set.
    typedef std::conditional_t<is_map_, std::pair<Kmer<k>, T_>, Kmer<k>> flat_type; // Underlying data type in the flat table.
    typedef std::conditional_t<is_map_, std::optional<T_>, bool> find_ret_t;    // Return type for find methods.
    static constexpr auto null_val = [](){ if constexpr(is_map_) return std::nullopt; else return false; }();   // Marker for absent keys (also values).

    std::size_t cap_;   // Capacity of the table in terms of buckets; adjusted to be a power of 2.

    flat_type* T;   // Flat table of size `cap x B`.

    // Metadata of the k-mers in a bucket.
    struct Metadata
    {
        uint8_t cs[B];  // Checksums of the k-mers.
        uint8_t min_coord[B];   // Minimizer coordinates within the k-mers.

        static constexpr uint8_t lock_mask = 0b1000'0000;   // Bit of `min_coord[0]` used to lock the bucket.
    };

    Metadata* M;    // Metadata table of the buckets.

    static constexpr uint64_t skew_th = B;
    std::conditional_t<is_map_, Hash_Map<Kmer<k>, T_>, Hash_Set<Kmer<k>>> overflow_;    // The overflow table in sequential setting.
    std::conditional_t<is_map_,
        boost::unordered::concurrent_flat_map<Kmer<k>, T_, Hash<Kmer<k>>>,
        boost::unordered::concurrent_flat_set<Kmer<k>, Hash<Kmer<k>>> > overflow_mt_;   // The overflow table in concurrent setting.

    uint64_t sum_len;   // Total sequence length provided.
    uint64_t sup_kmer_c;    // Number of super k-mers in the provided sequences to the table.

    // Following are for used for developer stats.
    ankerl::unordered_dense::map<uint64_t, ankerl::unordered_dense::set<Kmer<k>, Kmer_Hasher<k>>> min_to_kmer;  // Set of overflowing minimizers.


    // Returns `true` iff `key` matches to the key in `x`.
    static constexpr bool key_equals(const flat_type& x, const Kmer<k>& key) { if constexpr(is_map_) return x.first == key; else return x == key; }

    // Returns the size of the `b`'th bucket.
    std::size_t bucket_size(std::size_t b) const;

    // Locks the `b`'th bucket for exclusive-write access.
    void lock(std::size_t b);

    // Unlocks the `b`'th bucket from exclusive-write access.
    void unlock(std::size_t b);

    // Returns `true` iff the k-mer `key` with checksum `c` is present in the
    // `i`'th bucket.
    find_ret_t find(Kmer<k> key, uint8_t c, std::size_t i) const;

    // Directly places the k-mer `key` with checksum `c` into the `i`'th
    // bucket at index `j`.
    void place_at(const Kmer<k> key, const uint8_t c, const std::size_t i, const std::size_t j) { T[i * B + j] = key, M[i].cs[j] = c; }

    // Tries to insert the k-mer `key` with checksum `c` into the `i`'th bucket,
    // returning `true` iff `key` is present at the bucket after this operation.
    bool try_insert_at(Kmer<k> key, uint8_t c, std::size_t i);

    // Returns the overflow table. `mt_` denotes whether multiple threads can
    // access the table concurrently.
    constexpr auto& overflow() { if constexpr(mt_) return overflow_mt_; else return overflow_; }

public:

    // Constructs a streaming k-mer hash table to support up-to `max_sz` k-mers,
    // with a maximum load-factor of `lf`.
    Streaming_Kmer_Hash_Table(std::size_t max_sz, double lf = lf_default);

    // Deserializes a streaming k-mer hash table from the stream `is`.
    Streaming_Kmer_Hash_Table(std::ifstream& is);

    Streaming_Kmer_Hash_Table(const Streaming_Kmer_Hash_Table&) = delete;
    Streaming_Kmer_Hash_Table(Streaming_Kmer_Hash_Table&&) = delete;
    Streaming_Kmer_Hash_Table& operator=(const Streaming_Kmer_Hash_Table&) = delete;
    Streaming_Kmer_Hash_Table& operator=(Streaming_Kmer_Hash_Table&&) = delete;

    ~Streaming_Kmer_Hash_Table();

    // Returns the capacity (in terms of keys) of the table.
    auto capacity() const { return cap_ * B; }

    // Returns the size of the table.
    std::size_t size() const;

    // Returns the size of the overflow table.
    std::size_t overflow_size() const { return !mt_ ? overflow_.size() : overflow_mt_.size();  }

    // Clears the table.
    void clear();

    // Inserts the k-mer defined by the window `w` into the table.
    void insert(const Kmer_Window<k, l>& w);

    // Returns `true` iff the k-mer defined by the window `w` is present in the
    // table.
    find_ret_t find(const Kmer_Window<k, l>& w) const;

    // Serializes the table to the stream `os`.
    void serialize(std::ofstream& os) const;

    // Deserializes the table from the stream `is`.
    void deserialize(std::ifstream& is);

    // Reports developer statistics about the table.
    void report_stats() const;
};


// A window definining a k-mer with its associated metadata: its hash and
// minimizer.
template <uint16_t k, uint16_t l = 19>
class Kmer_Window
{
    template <uint16_t k_, bool mt_, typename T_, uint16_t l_> friend class Streaming_Kmer_Hash_Table;

    Directed_Vertex<k> v;
    Rolling_Hash<k, true> rh;
    Min_Iterator<k, l> min_it;

public:

    // Initializes the k-mer window at the beginning of the sequence `s`.
    void init(const char* const s)
    {
        v = Directed_Vertex<k>(Kmer<k>(s));
        rh.init(s);
        min_it.reset(s);
    }

    // Advances the window by one k-mer, by the character `ch`.
    void advance(const char ch)
    {
        assert(!DNA_Utility::is_placeholder(ch));

        v.roll_forward(DNA_Utility::map_base(ch));
        min_it.advance(ch);
        rh.advance(ch);
    }
};


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Streaming_Kmer_Hash_Table(const std::size_t max_sz, const double lf):
      cap_(ceil_pow_2((max_sz / lf) / B))
    , T(aligned_allocate<flat_type>(cap_ * B, L1_CACHE_LINE_SIZE))
    , M(aligned_allocate<Metadata>(cap_, L1_CACHE_LINE_SIZE))
{
    clear();
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Streaming_Kmer_Hash_Table(std::ifstream& is):
      cap_(0)
    , T(nullptr)
    , M(nullptr)
{
    deserialize(is);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline Streaming_Kmer_Hash_Table<k, mt_, T_, l>::~Streaming_Kmer_Hash_Table()
{
    deallocate(T), deallocate(M);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::clear()
{
    sum_len = 0, sup_kmer_c = 0;

    std::memset(M, 0, cap_ * sizeof(Metadata)); // Checksum `0` denotes absence of a key.

    overflow_.clear(), overflow_mt_.clear();

    min_to_kmer.clear();
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline std::size_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::bucket_size(const std::size_t b) const
{
    const auto zero  = _mm256_setzero_si256();
    const auto c_vec = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[b].cs));
    const auto cmp   = _mm256_cmpeq_epi8_mask(c_vec, zero);

    assert(std::countr_zero(cmp) + std::countl_one(cmp) == B);
    return std::countr_zero(cmp);   // This must hold if no deletion is performed at the table.
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::lock(const std::size_t b)
{
    while(true)
    {
        const auto exp_unlocked = M[b].min_coord[0] & ~Metadata::lock_mask;
        if(__sync_bool_compare_and_swap(&M[b].min_coord[0], exp_unlocked, exp_unlocked | Metadata::lock_mask))
            break;
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::unlock(const std::size_t b)
{
    assert(M[b].min_coord[0] & Metadata::lock_mask);
    M[b].min_coord[0] &= ~Metadata::lock_mask;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline std::size_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::size() const
{
    std::size_t sz = 0;
    for(std::size_t b = 0; b < cap_; ++b)
        sz += bucket_size(b);

    return sz + overflow_size();
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::find(const Kmer<k> key, const uint8_t c, const std::size_t i) const -> find_ret_t
{
    assert(i < cap_);

    const auto c_vec = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[i].cs));
    const auto query = _mm256_set1_epi8(static_cast<char>(c));
    auto cmp = _mm256_cmpeq_epi8_mask(c_vec, query);

    while(cmp)
    {
        const auto j = std::countr_zero(cmp);
        if(key_equals(T[i * B + j], key))
        {
            if constexpr(is_map_)   return T[i * B + j].second;
            else                    return true;
        }

        cmp &= (cmp - 1);
    }

    return null_val;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::try_insert_at(const Kmer<k> key, const uint8_t c, const std::size_t i)
{
    if constexpr(mt_)   lock(i);

    if(find(key, c, i))
    {
        if constexpr(mt_)   unlock(i);
        return true;
    }

    const auto c_vec = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[i].cs));
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


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::insert(const Kmer_Window<k, l>& w)
{
    const auto key = w.v.canonical();
    const auto h = w.min_it.hash();
    const auto idx_mask = cap_ - 1;

    const auto b = h & idx_mask;
    const auto c = std::max(w.rh.template checksum<8>(), uint64_t(1));  // Avoiding checksum 0 by overloading checksum 1.
    if(!try_insert_at(key, c, b))
    {
        // TODO: get the following hashes from the minimizer iterator (rolling hashes), instead of full-blown xxHash computations.
        const auto h_1 = XXH3_64bits_withSeed(&h, sizeof(h), 0);
        const auto h_2 = XXH3_64bits_withSeed(&h, sizeof(h), 1llu << 63);

        const auto s_1 = h_1 & idx_mask;
        const auto s_2 = h_2 & idx_mask;
        const auto p = std::min(s_1, s_2), q = std::max(s_1, s_2);
        assert(p != q);

        if constexpr(mt_)   lock(p), lock(q);
        if(find(key, c, p) || find(key, c, q))
        {
            if constexpr(mt_)   unlock(q), unlock(p);
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
                if constexpr(mt_)   unlock(q), unlock(p);
            }
            else if(sz_q < B)
            {
                place_at(key, c, q, sz_q);
                if constexpr(mt_)   unlock(q), unlock(p);
            }
            else
            {
                if constexpr(mt_)   unlock(q), unlock(p);
                overflow().insert(key);
            }
        }
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::find(const Kmer_Window<k, l>& w) const -> find_ret_t
{
    const auto key = w.v.canonical();
    const auto h = w.min_it.hash();
    const auto idx_mask = cap_ - 1;

    const auto b = h & idx_mask;
    const auto c = std::max(w.rh.template checksum<8>(), uint64_t(1));    // Avoiding checksum 0 by overloading checksum 1.
    if(const auto r = find(key, c, b))
        return r;

    // TODO: get the following hashes from the minimizer iterator (rolling hashes), instead of full-blown xxHash computations.
    const auto h_1 = XXH3_64bits_withSeed(&h, sizeof(h), 0);
    const auto h_2 = XXH3_64bits_withSeed(&h, sizeof(h), 1llu << 63);

    const auto p = h_1 & idx_mask;
    const auto q = h_2 & idx_mask;
    if(const auto r_p = find(key, c, p))    return r_p;
    if(const auto r_q = find(key, c, q))    return r_q;

    if(bucket_size(p) == B && bucket_size(q) == B)
        return overflow_.find(key);

    return null_val;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::serialize(std::ofstream& os) const
{
    ::serialize(k, os), ::serialize(l, os), ::serialize(B, os);

    ::serialize(cap_, os);
    ::serialize(T, cap_ * B, os), ::serialize(M, cap_, os);

    overflow_.serialize(os);

    ::serialize(sum_len, os), ::serialize(sup_kmer_c, os);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::deserialize(std::ifstream& is)
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

    if(T)
        deallocate(T), deallocate(M);

    T = aligned_allocate<Kmer<k>>(cap_ * B, L1_CACHE_LINE_SIZE), M = aligned_allocate<Metadata>(cap_, L1_CACHE_LINE_SIZE);
    ::deserialize(T, cap_ * B, is), ::deserialize(M, cap_, is);

    overflow_.deserialize(is);

    ::deserialize(sum_len, is), ::deserialize(sup_kmer_c, is);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::report_stats() const
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

    std::cerr   << "Complete table size: " << sz << ".\n";
    std::cerr   << "Overflow table size: " << overflow_.size() << " ("
                    << (overflow_.size() * 100.0 / sz) << "%).\n";
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
