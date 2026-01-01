
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
#include "boost/unordered/concurrent_flat_map.hpp"

#include <cstdint>
#include <cstddef>
#include <bit>
#include <atomic>
#include <optional>
#include <utility>
#include <type_traits>
#include <cstring>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <thread>
#include <algorithm>
#include <cassert>

#include <immintrin.h>


#define KACHE_HASH_LIKELY(cond)     __builtin_expect(cond, 1)
#define KACHE_HASH_UNLIKELY(cond)   __builtin_expect(cond, 0)


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
    typedef std::conditional_t<is_map_, std::pair<Kmer<k>, T_>, Kmer<k>> flat_t;    // Underlying data type in the flat table.
    typedef std::optional<T_> val_t;    // Type of values returned by hash-map operations.
    typedef std::conditional_t<is_map_, val_t, bool> find_ret_t;    // Return type for find / insert methods.
    static constexpr auto null_val = [](){ if constexpr(is_map_) return std::nullopt;   else return false; }(); // Marker for absent keys (also values).

    static constexpr auto& key(const flat_t& e) { if constexpr(is_map_) return e.first; else return e; }    // Returns the key from a table-element `e`.
    static constexpr auto& val(const flat_t& e) { if constexpr(is_map_) return e.second; else return e; }   // Returns the value from a table-element `e`.

    static constexpr auto min_orientation_mask = Min_Iterator<k, l>::min_orientation_mask();

    const double lf;    // Load factor.

    std::size_t cap_;   // Capacity of the table in terms of buckets; adjusted to be a power of 2.
    uint64_t resize_th; // At what size of the table to trigger the next resize.

    flat_t* T;  // Flat table of size `cap x B`.

    // Metadata of the k-mers in a bucket.
    struct Metadata
    {
        uint8_t cs[B];  // Checksums of the k-mers.
        uint8_t min_coord[B];   // Minimizer coordinates within the k-mers.

        static constexpr uint8_t lock_mask = 0b1000'0000;   // Bit of `min_coord[0]` used to lock the bucket.
    };

    Metadata* M;    // Metadata table of the buckets.

    flat_t* T_new;  // New flat table after a resize.
    Metadata* M_new;    // New metadata table of buckets after a resize.

    struct Overflow_Set_Val
    {
        uint8_t cs;
        uint8_t min_coord;

        Overflow_Set_Val(uint8_t cs, uint8_t min_coord): cs(cs), min_coord(min_coord) {}
    };

    struct Overflow_Map_Val
    {
        uint8_t cs;
        uint8_t min_coord;
        T_ val;

        Overflow_Map_Val(uint8_t cs, uint8_t min_coord, const T_& val): cs(cs), min_coord(min_coord), val(val) {}
    };

    static constexpr uint64_t skew_th = B;
    std::conditional_t<is_map_,
        boost::unordered::concurrent_flat_map<Kmer<k>, Overflow_Map_Val, Hash<Kmer<k>>>,
        boost::unordered::concurrent_flat_map<Kmer<k>, Overflow_Set_Val, Hash<Kmer<k>>> > overflow;    // The overflow table.

    std::atomic_uint16_t registered_thread_count;   // Number of threads registered to use the table.
    std::atomic_uint64_t approx_sz; // Approximate size (lower bound) of the hash table.
    std::vector<Padded<uint64_t>> checkpoint;   // `checkpoint[i]` is the number of elements remaining to be added by the `i`'th user thread before it checks for a resize.
    static constexpr uint64_t resize_checkpoint = 16384;    // Checkpoint number of elements for a user thread to add before checking for resize.

    uint64_t sum_len;   // Total sequence length provided.
    uint64_t sup_kmer_c;    // Number of super k-mers in the provided sequences to the table.

    // Following are for used for developer stats.
    ankerl::unordered_dense::map<uint64_t, ankerl::unordered_dense::set<Kmer<k>, Kmer_Hasher<k>>> min_to_kmer;  // Set of overflowing minimizers.


    // Returns `true` iff `key` matches to the key of the element `e`.
    static constexpr bool key_equals(const flat_t& e, const Kmer<k>& key) { return Streaming_Kmer_Hash_Table::key(e) == key; }

    // Returns the 32-bit equality mask of the vectors `x` and `y`.
    static uint32_t eq_mask(__m256i x, __m256i y);

    // Returns the size of the `b`'th bucket.
    std::size_t bucket_size(std::size_t b) const;

    // Locks the `b`'th bucket for exclusive-write access in concurrent
    // settings. `resize_` denotes whether the bucket requires to be locked for
    // resizing or not.
    template <bool resize_ = false> void lock(std::size_t b);

    // Locks the `p`'th and the `q`'th buckets for exclusive-write accesses in
    // concurrent settings, taking care of circular deadlocks and when the
    // buckets are the same. `resize_` denotes whether the buckets require to be
    // locked for resizing or not.
    template <bool resize_ = false> void lock_ordered(std::size_t p, std::size_t q);

    // Unlocks the `b`'th bucket from exclusive-write access in concurrent
    // settings. `resize_` denotes whether the bucket requires to be unlocked
    // during resizing or not.
    template <bool resize_ = false> void unlock(std::size_t b);

    // Unlocks the `p` and the `q`'th buckets from exclusive-write access in
    // concurrent settings, taking care of circular deadlocks and when the
    // buckets are the same. `resize_` denotes whether the buckets require to be
    // unlocked during resizing or not.
    template <bool resize_ = false> void unlock_ordered(std::size_t p, std::size_t q);

    // Returns the associated value iff the k-mer `key` with checksum `c` is
    // present in the `i`'th bucket.
    find_ret_t find(Kmer<k> key, uint8_t c, std::size_t i) const;

    // Returns the associated value iff the k-mer `key` is present in the
    // overflow table.
    find_ret_t find_in_overflow(Kmer<k> key) const;

    // Updates the associated value of the k-mer `key` with checksum `c` in the
    // `i`'th bucket through `f` iff it is present there, returning the earlier
    // value. Returns `null_opt` otherwise.
    template <typename F, typename V = T_>
    val_t find_and_update(Kmer<k> key, uint8_t c, std::size_t i, const F& f) requires (is_map_);

    // Directly places the element `e` with checksum `c` and minimizer-
    // coordinate `m` into the `i`'th bucket at index `j`.
    void place_at(const flat_t& e, const uint8_t c, const uint8_t m, const std::size_t i, const std::size_t j)
    { T[i * B + j] = e, M[i].cs[j] = c, M[i].min_coord[j] = m; }

    // Tries to insert the element `e` with checksum `c` and minimizer-
    // coordinate `m` into the `i`'th bucket. It succeeds iff the key was absent
    // in the bucket beforehand and there was empty space, returning `true` in
    // that case. Returns `false` otherwise; in this case if the key was already
    // in the bucket, it puts the associated value at `val`, or else puts
    // `null_val` at `val`.
    bool try_insert_at(const flat_t& e, uint8_t c, uint8_t m, std::size_t i, find_ret_t& val);

    // Tries to update the associated value of the k-mer `key` with checksum `c`
    // and minimizer-coordinate `m` at the `i`'th bucket through `f` iff it is
    // present there, and also puts the associated value to `old_val`. Otherwise
    // if the bucket has space, puts `key` with value `val` into it and puts
    // `null_val` at `old_val`. Returns `true` iff the upsert succeeds.
    template <typename F, typename V = T_>
    bool try_upsert_at(Kmer<k> key, uint8_t c, uint8_t m, std::size_t i, const F& f, const V& val, val_t& old_val) requires (is_map_);

    // Tries to insert the element `e` with checksum `c`, minimizer-hash `h`,
    // and -coordinate `m` into the table. It succeeds iff the key was absent in
    // the table beforehand, and returns `null_val` in that case. Otherwise it
    // returns the currently associated value.
    find_ret_t insert(const flat_t& e, uint8_t c, uint64_t h, uint8_t m);

    // Inserts the k-mer `key` with checksum `c`, minimizer-hash `h`, and
    // -coordinate `m` into the table with associated value `val` iff it is not
    // already present. Otherwise updates the currently associated value through
    // `f`. Returns the old associated value if the key was already present, and
    // returns `null_val` otherwise.
    template <typename F, typename V = T_>
    val_t upsert(const Kmer<k> key, uint8_t c, uint64_t h, uint8_t m, const F& f, const V& val) requires (is_map_);

    // Resizes the hash table to double its capacity.
    void resize();

public:

    class Token;

    // Constructs a streaming k-mer hash table to support up-to `max_sz` k-mers,
    // with a maximum load-factor of `lf`.
    Streaming_Kmer_Hash_Table(std::size_t max_sz = 1 << 22, double lf = lf_default);

    // Deserializes a streaming k-mer hash table from the stream `is`.
    // Streaming_Kmer_Hash_Table(std::ifstream& is);

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
    std::size_t overflow_size() const { return overflow.size();  }

    // Clears the table.
    void clear();

    // Inserts the k-mer defined by the window `w` into the table iff it is not
    // already present in the table, returning `false` in that case. Returns
    // `true` otherwise; i.e. returns the presence-status before the insertion.
    // The invoking user thread must provide its token `token`.
    bool insert(const Kmer_Window<k, l>& w, const Token& token) requires (!is_map_);

    // Inserts the k-mer defined by the window `w` and associated value `val`
    // into the table iff the key is not already present in the table, returning
    // `null_opt` in that case. Returns the associated current value otherwise;
    // i.e. returns the associated value pre-insert. The invoking user thread
    // must provide its token `token`.
    template <typename V = T_> val_t insert(const Kmer_Window<k, l>& w, const V& val, const Token& token) requires (is_map_);

    // Updates the associated value of the k-mer defined by `w` in the table
    // through `f`, iff it is present there. Otherwise inserts it into the table
    // with associated value `val`. Returns the old associated value if the key
    // was already present in the table, and returns `null_val` otherwise. The
    // invoking user thread must provide its token `token`.
    template <typename F, typename V = T_>
    val_t upsert(const Kmer_Window<k, l>& w, const F& f, const V& val, const Token& token) requires (is_map_);

    // Returns `true` iff the k-mer defined by the window `w` is present in the
    // table.
    find_ret_t find(const Kmer_Window<k, l>& w) const;

    // Registers a new user and returns a unique token for it.
    Token register_user() { return Token(registered_thread_count++); }

    class Iterator;

    // Returns an iterator to the beginning of the table.
    Iterator cbegin() const { return Iterator(*this, 0, 0); }

    // Returns an iterator to the (exclusive-)end of the table.
    Iterator cend() const { return Iterator(*this, cap_, -1); };

    // Serializes the table to the stream `os`.
    // void serialize(std::ofstream& os) const;

    // Deserializes the table from the stream `is`.
    // void deserialize(std::ifstream& is);

    // Reports developer statistics about the table.
    void report_stats() const;
};


// Token required to use the hash table.
template <uint16_t k, bool mt_, typename T_, uint16_t l>
class alignas(L1_CACHE_LINE_SIZE)
    Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Token
{
    friend class Streaming_Kmer_Hash_Table;

private:

    uint16_t id;    // Unique ID of the token.

    Token(const uint16_t id): id(id)
    {}

public:

    Token(): id(0)
    {}
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
class Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Iterator
{
    friend class Streaming_Kmer_Hash_Table<k, mt_, T_, l>;
    typedef Streaming_Kmer_Hash_Table<k, mt_, T_, l> ht_t;

private:

    const ht_t& ht; // Associated hash table.
    std::size_t b;  // Index of the bucket in the table where the iterator sits.
    int8_t slot;    // Index of the slot in the bucket where the iterator sits.

    // Constructs an iterator for the hash table `ht` that sits in bucket `b` at
    // slot `s`.
    Iterator(const ht_t& ht, const std::size_t b, const int8_t s): ht(ht), b(b), slot(s)
    {
        assert(b == 0 || b == ht.cap_);
        if(b < ht.cap_ && s >= 0 && s < static_cast<int8_t>(B) && !ht.M[b].cs[s])
            operator++();
    }

public:

    // Repositions the iterator to the next successive element in the map, or to
    // its end if there are no more elements.
    void operator++()
    {
        for(; b < ht.cap_; ++b)
        {
            for(slot++; slot < static_cast<int8_t>(B); ++slot)
                if(ht.M[b].cs[slot])
                    return;

            slot = -1;
        }
    }

    const ht_t::flat_t& operator*() { assert(b < ht.cap_ && 0 <= slot && slot < B); return ht.T[b * B + slot]; }

    auto checksum() const { assert(b < ht.cap_ && 0 <= slot && slot < B); return ht.M[b].cs[slot]; }

    auto min_coord() const { assert(b < ht.cap_ && 0 <= slot && slot < B); return ht.M[b].min_coord[slot]; }

    bool operator==(const Iterator& rhs) const { return &ht == &rhs.ht && b == rhs.b && slot == rhs.slot; }
};


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Streaming_Kmer_Hash_Table(const std::size_t max_sz, const double lf):
      lf(lf)
    , cap_(ceil_pow_2((max_sz / lf) / B))
    , resize_th(cap_ * lf)
    , T(aligned_allocate<flat_t>(cap_ * B, L1_CACHE_LINE_SIZE))
    , M(aligned_allocate<Metadata>(cap_, L1_CACHE_LINE_SIZE))
    , T_new(nullptr)
    , M_new(nullptr)
    , registered_thread_count(0)
    , approx_sz(0)
    , checkpoint(std::thread::hardware_concurrency(), resize_checkpoint)
{
    clear();
}


/*
template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Streaming_Kmer_Hash_Table(std::ifstream& is):
      cap_(0)
    , T(nullptr)
    , M(nullptr)
    , registered_thread_count(0)
{
    deserialize(is);
}
*/


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

    approx_sz = 0;
    std::for_each(checkpoint.begin(), checkpoint.end(), [](auto& c){ c = resize_checkpoint; });

    overflow.clear();

    min_to_kmer.clear();
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline uint32_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::eq_mask(const __m256i x, const __m256i y)
{
#ifdef AVX_512
    return _mm256_cmpeq_epi8_mask(x, y);
#else
    const auto cmp = _mm256_cmpeq_epi8(x, y);
    return static_cast<uint32_t>(_mm256_movemask_epi8(cmp));
#endif
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline std::size_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::bucket_size(const std::size_t b) const
{
    const auto zero  = _mm256_setzero_si256();
    const auto c_vec = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[b].cs));
    const auto cmp   = eq_mask(c_vec, zero);

    assert(std::countr_zero(cmp) + std::countl_one(cmp) == B);
    return std::countr_zero(cmp);   // This must hold if no deletion is performed at the table.
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool resize_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::lock(const std::size_t b)
{
    if constexpr(mt_)
        while(true)
        {
            if constexpr(!resize_)
            {
                const auto exp_unlocked = M[b].min_coord[0] & ~Metadata::lock_mask;
                if(__sync_bool_compare_and_swap(&M[b].min_coord[0], exp_unlocked, exp_unlocked | Metadata::lock_mask))
                    break;
            }
            else
            {
                const auto exp_unlocked = M_new[b].min_coord[0] & ~Metadata::lock_mask;
                if(__sync_bool_compare_and_swap(&M_new[b].min_coord[0], exp_unlocked, exp_unlocked | Metadata::lock_mask))
                    break;
            }
        }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool resize_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::lock_ordered(const std::size_t p, const std::size_t q)
{
    if constexpr(mt_)
        KACHE_HASH_LIKELY(p != q) ?
            lock<resize_>(std::min(p, q)), lock<resize_>(std::max(p, q)) :
            lock<resize_>(p);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool resize_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::unlock(const std::size_t b)
{
    if constexpr(mt_)
    {
        if constexpr(!resize_)
        {
            assert(M[b].min_coord[0] & Metadata::lock_mask);
            M[b].min_coord[0] &= ~Metadata::lock_mask;
        }
        else
        {
            assert(M_new[b].min_coord[0] & Metadata::lock_mask);
            M_new[b].min_coord[0] &= ~Metadata::lock_mask;
        }
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool resize_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::unlock_ordered(const std::size_t p, const std::size_t q)
{
    if constexpr(mt_)
        KACHE_HASH_LIKELY(p != q) ?
                unlock<resize_>(std::max(p, q)), unlock<resize_>(std::min(p, q)) :
                unlock<resize_>(p);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline std::size_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::size() const
{
    std::size_t sz = approx_sz;
    std::for_each(checkpoint.cbegin(), checkpoint.cend(), [&](auto& c){ sz += c.unwrap() > 0 ? resize_checkpoint - c.unwrap() : 0; });
    return sz;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::find(const Kmer<k> key, const uint8_t c, const std::size_t i) const -> find_ret_t
{
    assert(i < cap_);

    const auto c_vec = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[i].cs));
    const auto query = _mm256_set1_epi8(static_cast<char>(c));
    const auto cmp = eq_mask(c_vec, query);

    const auto j = cmp ? __builtin_ctz(cmp) : 0;
    if(M[i].cs[j] && key_equals(T[i * B + j], key))
    {
        if constexpr(is_map_)   return T[i * B + j].second;
        else                    return true;
    }

    auto cmp_rem = cmp & (cmp - 1);
    while(cmp_rem)
    {
        const auto j = __builtin_ctz(cmp_rem);
        if(key_equals(T[i * B + j], key))
        {
            if constexpr(is_map_)   return T[i * B + j].second;
            else                    return true;
        }

        cmp_rem &= (cmp_rem - 1);
    }

    return null_val;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::find_in_overflow(const Kmer<k> key) const -> find_ret_t
{
    if constexpr(!is_map_)
        return overflow.count(key);
    else
    {
        auto val = null_val;
        overflow.cvisit(key, [&](const auto& kv){ val = kv.second.val; });
        return val;
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename F, typename V>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::find_and_update(const Kmer<k> key, const uint8_t c, const std::size_t i, const F& f) -> val_t requires (is_map_)
{
    assert(i < cap_);

    const auto c_vec = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[i].cs));
    const auto query = _mm256_set1_epi8(static_cast<char>(c));
    auto cmp = eq_mask(c_vec, query);

    while(cmp)
    {
        const auto j = __builtin_ctz(cmp);
        if(key_equals(T[i * B + j], key))
        {
            const auto prev_val = T[i * B + j].second;
            T[i * B + j].second = f(prev_val);
            return prev_val;
        }

        cmp &= (cmp - 1);
    }

    return null_val;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::try_insert_at(const flat_t& e, const uint8_t c, const uint8_t m, const std::size_t i, find_ret_t& val)
{
    lock(i);

    const auto& key = this->key(e);
    if(const auto r = find(key, c, i))
    {
        unlock(i);
        val = r;
        return false;
    }

    const auto c_vec = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[i].cs));
    const auto zero = _mm256_setzero_si256();
    const auto cmp = eq_mask(c_vec, zero);

    const auto j = std::countr_zero(cmp);
    if(j == B)
    {
        unlock(i);
        val = null_val;
        return false;
    }

    place_at(e, c, m, i, j);
    unlock(i);
    return true;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename F, typename V>
inline bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::try_upsert_at(const Kmer<k> key, const uint8_t c, const uint8_t m, const std::size_t i, const F& f, const V& val, val_t& old_val) requires (is_map_)
{
    lock(i);

    if(const auto r = find_and_update(key, c, i, f))
    {
        unlock(i);
        old_val = r;
        return true;
    }

    const auto c_vec = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[i].cs));
    const auto zero = _mm256_setzero_si256();
    const auto cmp = eq_mask(c_vec, zero);

    const auto j = std::countr_zero(cmp);
    if(j == B)
    {
        unlock(i);
        return false;
    }

    place_at(std::make_pair(key, val), c, m, i, j);
    unlock(i);

    old_val = null_val;
    return true;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::insert(const Kmer_Window<k, l>& w, const Token& token) requires (!is_map_)
{
    uint8_t m;
    const auto c = std::max(w.rh.template checksum<8>(), uint64_t(1));  // Avoiding checksum 0 by overloading checksum 1.
    const auto h = w.min_it.hash(m);
    m = w.v.in_canonical_form() ? m : (m ^ min_orientation_mask);
    const auto r = insert(w.v.canonical(), c, h, m);

    auto& cp = checkpoint[token.id].unwrap();
    cp -= (r == null_val);
    if(KACHE_HASH_UNLIKELY(cp == 0))
    {
        approx_sz += resize_checkpoint;
        cp = resize_checkpoint;
        // TODO: check resize.
    }

    return r;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename V>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::insert(const Kmer_Window<k, l>& w, const V& val, const Token& token) -> val_t requires (is_map_)
{
    uint8_t m;
    const auto c = std::max(w.rh.template checksum<8>(), uint64_t(1));  // Avoiding checksum 0 by overloading checksum 1.
    const auto h = w.min_it.hash(m);
    m = w.v.in_canonical_form() ? m : (m ^ min_orientation_mask);
    const auto r = insert(std::make_pair(w.v.canonical(), val), c, h, m);

    auto& cp = checkpoint[token.id].unwrap();
    cp -= (r == null_val);
    if(KACHE_HASH_UNLIKELY(cp == 0))
    {
        approx_sz += resize_checkpoint;
        cp = resize_checkpoint;
        // TODO: check resize.
    }

    return r;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename F, typename V>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::upsert(const Kmer_Window<k, l>& w, const F& f, const V& val, const Token& token) -> val_t requires (is_map_)
{
    (void)token;
    uint8_t m;
    const auto c = std::max(w.rh.template checksum<8>(), uint64_t(1));  // Avoiding checksum 0 by overloading checksum 1.
    const auto h = w.min_it.hash(m);
    m = w.v.in_canonical_form() ? m : (m ^ min_orientation_mask);
    const auto r = upsert(w.v.canonical(), c, h, m, f, val);

    auto& cp = checkpoint[token.id].unwrap();
    cp -= (r == null_val);
    if(KACHE_HASH_UNLIKELY(cp == 0))
    {
        approx_sz += resize_checkpoint;
        cp = resize_checkpoint;
        // TODO: check resize.
    }

    return r;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::insert(const flat_t& x, const uint8_t c, const uint64_t h, const uint8_t m) -> find_ret_t
{
    const auto idx_mask = cap_ - 1;

    const auto b = h & idx_mask;
    find_ret_t v;
    if(const auto r = try_insert_at(x, c, m, b, v))
        return null_val;
    else if(v)
        return v;
    else
    {
        // TODO: get the following hashes from the minimizer iterator (rolling hashes), instead of full-blown xxHash computations.
        const auto h_1 = XXH3_64bits_withSeed(&h, sizeof(h), 0);
        const auto h_2 = XXH3_64bits_withSeed(&h, sizeof(h), 1llu << 63);

        const auto s_1 = h_1 & idx_mask;
        const auto s_2 = h_2 & idx_mask;
        const auto p = std::min(s_1, s_2), q = std::max(s_1, s_2);
        assert(p != q);

        const auto key = this->key(x);

        lock_ordered(p, q);
        if(const auto r_p = find(key, c, p))
        {
            unlock_ordered(p, q);
            return r_p;
        }
        else if(const auto r_q = find(key, c, q))
        {
            unlock_ordered(p, q);
            return r_q;
        }
        else
        {
            const auto zero = _mm256_setzero_si256();

            const auto c_vec_p = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[p].cs));
            const auto c_vec_q = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[q].cs));

            const auto cmp_p = eq_mask(c_vec_p, zero);
            const auto cmp_q = eq_mask(c_vec_q, zero);

            const auto sz_p = static_cast<std::size_t>(std::countr_zero(cmp_p));
            const auto sz_q = static_cast<std::size_t>(std::countr_zero(cmp_q));
            assert(std::countr_zero(cmp_p) + std::countl_one(cmp_p) == B), assert(std::countr_zero(cmp_q) + std::countl_one(cmp_q) == B);

            if(sz_p < sz_q)
            {
                place_at(x, c, m, p, sz_p);
                unlock_ordered(p, q);
            }
            else if(sz_q < B)
            {
                place_at(x, c, m, q, sz_q);
                unlock_ordered(p, q);
            }
            else
            {
                unlock_ordered(p, q);
                if constexpr(!is_map_)  // TODO: incorrect value for now, due to limitations in the 3rd-party map APIs.
                    return overflow.insert(std::pair(x, Overflow_Set_Val(c, m))) ? null_val : find_ret_t(true);
                else
                    return overflow.insert(std::pair(x.first, Overflow_Map_Val(c, m, x.second))) ? null_val : find_ret_t(true);
            }

            return null_val;
        }
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename F, typename V>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::upsert(const Kmer<k> key, const uint8_t c, const uint64_t h, const uint8_t m, const F& f, const V& val) -> val_t requires (is_map_)
{
    const auto idx_mask = cap_ - 1;

    const auto b = h & idx_mask;
    val_t old_val;
    if(try_upsert_at(key, c, m, b, f, val, old_val))
        return old_val;
    else
    {
        // TODO: get the following hashes from the minimizer iterator (rolling hashes), instead of full-blown xxHash computations.
        const auto h_1 = XXH3_64bits_withSeed(&h, sizeof(h), 0);
        const auto h_2 = XXH3_64bits_withSeed(&h, sizeof(h), 1llu << 63);

        const auto p = h_1 & idx_mask;
        const auto q = h_2 & idx_mask;

        lock_ordered(p, q);
        if(const auto r_p = find_and_update(key, c, p, f))
        {
            unlock_ordered(p, q);
            return r_p;
        }
        else if(const auto r_q = find_and_update(key, c, q, f))
        {
            unlock_ordered(p, q);
            return r_q;
        }
        else
        {
            const auto zero = _mm256_setzero_si256();

            const auto c_vec_p = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[p].cs));
            const auto c_vec_q = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[q].cs));

            const auto cmp_p = eq_mask(c_vec_p, zero);
            const auto cmp_q = eq_mask(c_vec_q, zero);

            const auto sz_p = static_cast<std::size_t>(std::countr_zero(cmp_p));
            const auto sz_q = static_cast<std::size_t>(std::countr_zero(cmp_q));
            assert(std::countr_zero(cmp_p) + std::countl_one(cmp_p) == B), assert(std::countr_zero(cmp_q) + std::countl_one(cmp_q) == B);

            if(sz_p < sz_q)
            {
                place_at(std::make_pair(key, val), c, m, p, sz_p);
                unlock_ordered(p, q);
            }
            else if(sz_q < B)
            {
                place_at(std::make_pair(key, val), c, m, q, sz_q);
                unlock_ordered(p, q);
            }
            else
            {
                unlock_ordered(p, q);

                val_t old_val = null_val;
                overflow.emplace_or_visit(key, Overflow_Map_Val(c, m, val), [&](auto& kv)
                {
                    old_val = kv.second.val;
                    kv.second.val = f(kv.second.val);
                });

                return old_val;
            }

            return null_val;
        }
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::resize()
{
    const auto old_cap = cap_;
    cap_ *= 2;
    T_new = aligned_allocate<flat_t>(cap_ * B, L1_CACHE_LINE_SIZE);
    M_new = aligned_allocate<Metadata>(cap_, L1_CACHE_LINE_SIZE);
    std::memset(M_new, 0, cap_ * sizeof(Metadata)); // Checksum `0` denotes absence of a key.

    for(std::size_t b = 0; b < old_cap; ++b)
        for(std::size_t slot = 0; slot < B; ++slot)
            if(M[b].cs[slot])
            {
                const auto kv = T[b * B + slot];
                const auto key = this->key(kv);
                const auto c = M[b].cs[slot];
                const auto m = M[b].min_coord[slot];

                const uint8_t min_idx = m & Min_Iterator<k, l>::min_index_mask();
                const bool in_fwd = m & Min_Iterator<k, l>::min_orientation_mask();
                const auto min = in_fwd ?
                                    key.template lmer_at<l>(min_idx) : Kmer_Utility::reverse_complement<l>(key.template lmer_at<l>(k - min_idx - l));
                const auto h = Min_Iterator<k, l>::lmer_hash(min);

                // TODO: insert item into the new table.
            }

    // TODO: insert items from the overflow table into the new table.
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

    if(bucket_size(b) < B)
        return null_val;

    // TODO: get the following hashes from the minimizer iterator (rolling hashes), instead of full-blown xxHash computations.
    const auto h_1 = XXH3_64bits_withSeed(&h, sizeof(h), 0);
    const auto h_2 = XXH3_64bits_withSeed(&h, sizeof(h), 1llu << 63);

    const auto p = h_1 & idx_mask;
    const auto q = h_2 & idx_mask;
    if(const auto r_p = find(key, c, p))    return r_p;
    if(const auto r_q = find(key, c, q))    return r_q;

    if(bucket_size(p) == B && bucket_size(q) == B)
        return find_in_overflow(key);

    return null_val;
}


/*
template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::serialize(std::ofstream& os) const
{
    ::serialize(k, os), ::serialize(l, os), ::serialize(B, os);

    ::serialize(cap_, os);
    ::serialize(T, cap_ * B, os), ::serialize(M, cap_, os);

    // overflow_.serialize(os);

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

    // overflow_.deserialize(is);

    ::deserialize(sum_len, is), ::deserialize(sup_kmer_c, is);
}
*/


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
    std::cerr   << "Overflow table size: " << overflow.size() << " ("
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
