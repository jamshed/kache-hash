
#ifndef STREAMING_KMER_HASHSET_HPP
#define STREAMING_KMER_HASHSET_HPP



#include "DNA_Utility.hpp"
#include "Kmer.hpp"
#include "Kmer_Utility.hpp"
#include "Kmer_Batch.hpp"
#include "Concurrent_Hash_Table.hpp"
#include "RW_Lock.hpp"
#include "utility.hpp"

#include <cstdint>
#include <cstddef>
#include <bit>
#include <atomic>
#include <tuple>
#include <optional>
#include <utility>
#include <type_traits>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <thread>
#include <algorithm>
#include <cassert>

#include <immintrin.h>


namespace kache_hash
{


#define KACHE_HASH_LIKELY(cond)     __builtin_expect(cond, 1)
#define KACHE_HASH_UNLIKELY(cond)   __builtin_expect(cond, 0)


// Hash table particularly suited for streaming `k`-mer operations, for
// `k <= 32`. `mt_` denotes whether a concurrent table is desired. If a hash map
// is required, the mapped type is `T_`. `l`-minimizers are used in underlying
// computations.
// NB: concurrent insertions are guaranteed to function correctly only if
// deletions are not performed, as of now.
template <uint16_t k, bool mt_, typename T_ = void, uint16_t l = 19>
class Streaming_Kmer_Hash_Table
{
    static_assert(k <= 32, "The k-mers need to fit in a machine word.");

public:

    class Token;

    static constexpr uint64_t max_user = 128;   // Maximum number of users allowed to use the table: a thread per core of a 2×64-core node, for now.

private:

    static constexpr double lf_default_ = 0.8;  // Default load factor: 80%.
    static constexpr double of_default = 0.1;   // Fixed overflow factor: 10%.
    static constexpr std::size_t B = 32;    // Bucket size: a bucket's checksums fill a 256-bit vector.

    typedef uint32_t mask_t;    // A bit per slot of a bucket.

    static constexpr bool is_set_ = std::is_void_v<T_>; // Whether a hash-set.
    static constexpr bool is_map_ = !is_set_;   // Whether a hash-map.
    typedef std::conditional_t<is_map_, std::pair<Kmer<k>, T_>, Kmer<k>> flat_t;    // Underlying data type in the flat table.
    typedef std::optional<T_> val_t;    // Type of values returned by hash-map operations.
    typedef std::conditional_t<is_map_, val_t, bool> find_ret_t;    // Return type for find methods.
    static constexpr auto null_val = [](){ if constexpr(is_map_) return std::nullopt;   else return false; }(); // Marker for absent keys (also values).

    static constexpr auto& key(const flat_t& e) { if constexpr(is_map_) return e.first; else return e; }    // Returns the key from a table-element `e`.

    const double lf;    // Load factor.

    std::size_t cap_;   // Capacity of the table in terms of buckets; adjusted to be a power of 2.
    std::atomic_uint64_t resize_th; // At what size of the table to trigger the next resize.

    // Metadata of the k-mers in a bucket.
    struct Metadata
    {
        uint8_t cs[B];  // Checksums of the k-mers.
        uint8_t min_coord[B];   // Minimizer coordinates within the k-mers.

        static constexpr uint8_t lock_mask = 0b1000'0000;   // Bit of `min_coord[0]` used to lock the bucket.
    };

    struct Overflow_Set_Val
    {
        uint8_t cs;
        uint8_t min_coord;

        Overflow_Set_Val() {}
        Overflow_Set_Val(uint8_t cs, uint8_t min_coord): cs(cs), min_coord(min_coord) {}
    };

    struct Overflow_Map_Val
    {
        uint8_t cs;
        uint8_t min_coord;
        T_ val;

        Overflow_Map_Val() {}
        Overflow_Map_Val(uint8_t cs, uint8_t min_coord, const T_& val): cs(cs), min_coord(min_coord), val(val) {}
    };

    typedef std::conditional_t<is_set_,
        Concurrent_Hash_Table<Kmer<k>, Overflow_Set_Val, Hash<Kmer<k>>>,
        Concurrent_Hash_Table<Kmer<k>, Overflow_Map_Val, Hash<Kmer<k>>> > ov_t; // The overflow table type.


    flat_t* T;  // Flat table of size `cap x B`.
    Metadata* M;    // Metadata table of the buckets.
    ov_t* ov;   // The overflow table.

    flat_t* T_new;  // New flat table after a resize.
    Metadata* M_new;    // New metadata table of buckets after a resize.
    ov_t* ov_new;    // New overflow table after a resize.

    // Returns the metadata of the `i`'th bucket.
    Metadata& md(const std::size_t i) const { return M[i]; }

    // Returns the `j`'th slot of the `i`'th bucket.
    flat_t& slot(const std::size_t i, const std::size_t j) const { return T[i * B + j]; }

    // Returns the metadata of the `i`'th new bucket in a resize.
    Metadata& md_new(const std::size_t i) const { return M_new[i]; }

    // Returns the `j`'th slot of the `i`'th new bucket in a resize.
    flat_t& slot_new(const std::size_t i, const std::size_t j) const { return T_new[i * B + j]; }

    std::atomic_uint16_t registered_thread_count;   // Number of threads registered to use the table.

    std::atomic_uint64_t approx_sz; // Size of the hash table, short of the user-local counts not yet flushed into it.
    std::vector<Padded<uint64_t>> checkpoint;   // `checkpoint[i]` is the number of keys the `i`'th user has added and not yet flushed into `approx_sz`.
    static constexpr uint64_t resize_checkpoint = 16384;    // Most keys a user adds before flushing its count and checking for a resize.
    std::atomic_uint64_t checkpoint_step;   // Number of keys a user adds before flushing its count and checking for a resize: unflushed keys stay under an eighth of the resize threshold.
    RW_Lock<max_user> table_lock; // Readers-writers lock to resize the table exclusively.

    const uint64_t resize_worker_c; // Number of worker threads to spawn for a concurrent resize.


    // Returns `true` iff `key` matches to the key of the element `e`.
    static constexpr bool key_equals(const flat_t& e, const Kmer<k>& key) { return Streaming_Kmer_Hash_Table::key(e) == key; }

    typedef __m256i cs_vec; // A vector of the checksums of a bucket.

    // Returns the checksums of the bucket with metadata `md`.
    static cs_vec load_cs(const auto& md);

    // Returns the checksum `c` broadcast to a vector.
    static cs_vec splat_cs(uint8_t c);

    // Returns the vector of empty checksums.
    static cs_vec zero_cs() { return splat_cs(0); }

    // Returns the mask of the slots whose checksums are equal in `x` and `y`.
    static mask_t eq_mask(cs_vec x, cs_vec y);

    // Returns the checksums `v` with the `j`'th one set to `c`, which is empty
    // in `v`.
    static cs_vec set_cs(cs_vec v, std::size_t j, uint8_t c);

    // Returns the size of the `b`'th bucket.
    std::size_t bucket_size(std::size_t b) const;

    // Returns the pair of the numbers of keys at bucket `b` which have `b` as
    // their primary bucket and as their secondary bucket, respectively.
    std::pair<std::size_t, std::size_t> distribution(std::size_t b) const;

    // Locks the `b`'th bucket for exclusive-write access in concurrent
    // settings, i.e. in a concurrent table, or in a resize with multiple
    // workers. `resize_` denotes whether the bucket requires to be locked for
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

    // Locks the distinct ones of the buckets `x`, `y`, and `z` in increasing
    // order of index, in concurrent settings.
    void lock_sorted(std::size_t x, std::size_t y, std::size_t z);

    // Unlocks the distinct ones of the buckets `x`, `y`, and `z`, in
    // concurrent settings.
    void unlock_sorted(std::size_t x, std::size_t y, std::size_t z);

    // Returns the associated value iff the k-mer `key` is present in the
    // overflow table.
    __attribute__((noinline)) find_ret_t find_in_overflow(Kmer<k> key) const;

    // Directly places the element `e` with checksum `c` and minimizer-
    // coordinate `m` into the `i`'th bucket at index `j`. The bucket's lock-bit
    // in `min_coord[0]` is kept as is.
    void place_at(const flat_t& e, uint8_t c, uint8_t m, std::size_t i, std::size_t j);

    // Directly places the element `e` with checksum `c` and minimizer-
    // coordinate `m` into the `i`'th bucket at index `j` of the new resized
    // table. The bucket's lock-bit in `min_coord[0]` is kept as is.
    void place_at_resize(const flat_t& e, uint8_t c, uint8_t m, std::size_t i, std::size_t j);

    // Tries to insert the element `x` with checksum `c` and minimizer-
    // coordinate `m` into the `i`'th bucket of the newly resized table. It
    // succeeds iff there was empty space, returning `true` in that case.
    // Returns `false` otherwise.
    bool try_insert_at_resize(const flat_t& x, uint8_t c, uint8_t m, std::size_t i);

    // Inserts the element `x` with checksum `c` and minimizer-coordinate `m`
    // into the resized table.
    void insert_at_resize(const flat_t& x, uint8_t c, uint8_t m);

    // Returns the minimizer hash of the canonical k-mer `key` with the
    // minimizer-coordinate `m`.
    static uint64_t minimizer_hash(Kmer<k> key, uint8_t m);

    // Returns the secondary buckets of the minimizer hash `h`.
    __attribute__((noinline)) std::pair<std::size_t, std::size_t> secondaries(uint64_t h) const;

    // Returns the number of keys a user adds before flushing its count and
    // checking for a resize, for the resize threshold `th`.
    static uint64_t step_for(uint64_t th) { return std::clamp<uint64_t>(th / (8 * max_user), 1, resize_checkpoint); }

    // Counts a key newly added by the `id`'th user. Returns whether the table
    // is due for a resize.
    bool add_key(uint16_t id);

    // Tries to resize the hash table to double its capacity: does so iff it is
    // due.
    __attribute__((noinline, cold)) void try_resize();

    // Resizes the hash table to double its capacity.
    void resize();

    // Descriptor of a k-mer for the streaming operations: what the table
    // operations need of it.
    struct Kmer_Desc
    {
        Kmer<k> key;    // Canonical form of the k-mer.
        uint64_t h; // Hash of its minimizer.
        uint8_t c;  // Its checksum.
        uint8_t m;  // Coordinate of its minimizer, oriented per the canonical form.
    };

    // Applies `op(d, q, i)` to the descriptor `d` of each k-mer of the
    // sequences `seqs[0..c)` of lengths `lens[0..c)`, the k-mer being at index
    // `i` of the `q`'th sequence, in order. Non-ACGT characters split the
    // sequences. The k-mers of a stretch of ACGT bases are described in chunks.
    template <typename Op_>
    static void stream(const char* const* seqs, const std::size_t* lens, std::size_t c, const Op_& op);

    // Returns the index of the slot holding the k-mer `key` with checksum `c`
    // in the `i`'th bucket, or `B` if it is absent there. Puts the index of the
    // bucket's first empty slot at `e`, `B` if it is full.
    uint32_t probe(Kmer<k> key, uint8_t c, std::size_t i, uint32_t& e) const { return probe(key, c, i, load_cs(md(i)), e); }

    // Returns the index of the slot holding the k-mer `key` with checksum `c`
    // in the `i`'th bucket, which has the checksums `c_vec`; or `B` if it is
    // absent there. Puts the index of the bucket's first empty slot at `e`, `B`
    // if it is full.
    uint32_t probe(Kmer<k> key, uint8_t c, std::size_t i, cs_vec c_vec, uint32_t& e) const;

    // Inserts the k-mer `key` with checksum `c` and coordinate `m` into the
    // overflow table iff it is absent there, with value `val` in a map; also in
    // a map, updates the value `v` of a present one to `f(v)`. Returns whether
    // it was absent.
    // NB: Out of line, as the streaming updates rarely reach the overflow
    // table.
    template <typename F, typename V>
    __attribute__((noinline)) bool ov_update(Kmer<k> key, uint8_t c, uint8_t m, const F& f, const V& val);

    // Updates the table with the k-mers of the sequences `seqs[0..c)` of
    // lengths `lens[0..c)`: inserts each absent k-mer with value `val`, and in
    // a map, updates the value `v` of each present one to `f(v)`. Returns the
    // number of k-mers inserted. The invoking user thread must provide its
    // token `token`.
    // NB: The streaming operations force-inline all their callees but the cold
    // ones (resize), rather than leaving it to the compiler's heuristics, whose
    // choices flip with small changes in these large functions.
    template <typename F, typename V>
    std::size_t update(const char* const* seqs, const std::size_t* lens, std::size_t c, const F& f, const V& val, const Token& token);

public:

    // Constructs a streaming k-mer hash table to support up-to `max_sz` k-mers,
    // with a maximum load-factor of `lf`. Resizes use `resize_worker_c`
    // threads.
    Streaming_Kmer_Hash_Table(std::size_t max_sz = 1 << 22, uint64_t resize_worker_c = 16, double lf = lf_default_);

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
    std::size_t overflow_size() const { return ov->size();  }

    // Inserts the k-mers of the sequence `seq` of length `len` into the set,
    // skipping the ones already present. Non-ACGT characters split the
    // sequence. Returns the number of k-mers inserted. The invoking user
    // thread must provide its token `token`.
    std::size_t insert(const char* seq, std::size_t len, const Token& token) requires (is_set_) { return update(&seq, &len, 1, [](){}, 0, token); }

    // Inserts the k-mers of the sequences `seqs[0..c)` of lengths `lens[0..c)`
    // into the set, as `insert` does for one sequence. Batching sequences,
    // such as reads, lets their k-mers be described at once.
    std::size_t insert(const char* const* seqs, const std::size_t* lens, std::size_t c, const Token& token) requires (is_set_) { return update(seqs, lens, c, [](){}, 0, token); }

    // Updates the associated value `v` of each k-mer of the sequence `seq` of
    // length `len` present in the map to `f(v)`, and inserts each absent one
    // with associated value `val`, one k-mer at a time, in an unspecified
    // order. Non-ACGT characters split the sequence. Returns the number of
    // k-mers inserted. The invoking user thread must provide its token
    // `token`.
    template <typename F, typename V = T_> std::size_t upsert(const char* seq, std::size_t len, const F& f, const V& val, const Token& token) requires (is_map_) { return update(&seq, &len, 1, f, val, token); }

    // Upserts the k-mers of the sequences `seqs[0..c)` of lengths `lens[0..c)`
    // into the map, as `upsert` does for one sequence. Batching sequences,
    // such as reads, lets their k-mers be described at once.
    template <typename F, typename V = T_> std::size_t upsert(const char* const* seqs, const std::size_t* lens, std::size_t c, const F& f, const V& val, const Token& token) requires (is_map_) { return update(seqs, lens, c, f, val, token); }

    // Looks up the k-mers of the sequence `seq` of length `len`, calling
    // `g(i, r)` for the k-mer at index `i`, with `r` as `true` in a set, or
    // its value in a map, iff it is present, in an unspecified order.
    // Non-ACGT characters split the sequence, and their k-mers are skipped.
    // Must not overlap updates that can resize the table.
    template <typename G> void find(const char* seq, std::size_t len, const G& g) const { find(&seq, &len, 1, [&](std::size_t, const std::size_t i, const find_ret_t& r){ g(i, r); }); }

    // Looks up the k-mers of the sequences `seqs[0..c)` of lengths
    // `lens[0..c)`, calling `g(q, i, r)` for the k-mer at index `i` of the
    // `q`'th sequence, as `find` does for one sequence.
    template <typename G> void find(const char* const* seqs, const std::size_t* lens, std::size_t c, const G& g) const;

    // Returns the number of k-mers of the sequence `seq` of length `len`
    // present in the table. Must not overlap updates that can resize the
    // table.
    std::size_t count(const char* seq, std::size_t len) const { return count(&seq, &len, 1); }

    // Returns the number of k-mers of the sequences `seqs[0..c)` of lengths
    // `lens[0..c)` present in the table. Must not overlap updates that can
    // resize the table.
    std::size_t count(const char* const* seqs, const std::size_t* lens, std::size_t c) const;

    // Registers a new user and returns a unique token for it.
    Token register_user() { return Token(registered_thread_count++); }

    // Returns the percentages of keys in the table that are placed into their
    // primary bucket, their secondary buckets, and the overflow table.
    std::tuple<double, double, double> key_distribution() const;
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


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Streaming_Kmer_Hash_Table(const std::size_t max_sz, const uint64_t resize_worker_c, const double lf):
      lf(lf)
    , cap_(ceil_pow_2(std::max<std::size_t>(1, (max_sz / lf) / B)))
    , resize_th(capacity() * lf)
    , T(aligned_allocate<flat_t>(cap_ * B, L1_CACHE_LINE_SIZE))
    , M(aligned_allocate<Metadata>(cap_, L1_CACHE_LINE_SIZE))
    , ov(new ov_t(capacity() * of_default))
    , T_new(nullptr)
    , M_new(nullptr)
    , ov_new(nullptr)
    , registered_thread_count(0)
    , approx_sz(0)
    , checkpoint(max_user, Padded<uint64_t>(0))
    , checkpoint_step(step_for(resize_th))
    , resize_worker_c(resize_worker_c)
{
    std::memset(M, 0, cap_ * sizeof(Metadata)); // Checksum `0` denotes absence of a key.
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline Streaming_Kmer_Hash_Table<k, mt_, T_, l>::~Streaming_Kmer_Hash_Table()
{
    deallocate(T), deallocate(M), delete(ov);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::load_cs(const auto& md) -> cs_vec
{
    return _mm256_load_si256(reinterpret_cast<const __m256i*>(md.cs));
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::splat_cs(const uint8_t c) -> cs_vec
{
    return _mm256_set1_epi8(static_cast<char>(c));
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::eq_mask(const cs_vec x, const cs_vec y) -> mask_t
{
    return static_cast<uint32_t>(_mm256_movemask_epi8(_mm256_cmpeq_epi8(x, y)));
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::set_cs(const cs_vec v, const std::size_t j, const uint8_t c) -> cs_vec
{
    const auto iota = _mm256_setr_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31);
    return _mm256_or_si256(v, _mm256_and_si256(_mm256_cmpeq_epi8(iota, _mm256_set1_epi8(static_cast<char>(j))), _mm256_set1_epi8(static_cast<char>(c))));
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline std::size_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::bucket_size(const std::size_t b) const
{
    const auto zero  = zero_cs();
    const auto c_vec = load_cs(md(b));
    const auto cmp   = eq_mask(c_vec, zero);

    assert(std::countr_zero(cmp) + std::countl_one(cmp) == B);
    return std::countr_zero(cmp);   // This must hold if no deletion is performed at the table.
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline std::pair<std::size_t, std::size_t> Streaming_Kmer_Hash_Table<k, mt_, T_, l>::distribution(const std::size_t b) const
{
    const auto idx_mask = cap_ - 1;

    std::size_t p = 0, s = 0;
    for(std::size_t j = 0; j < B; ++j)
        if(md(b).cs[j])
        {
            const auto h = minimizer_hash(key(slot(b, j)), md(b).min_coord[j] & ~Metadata::lock_mask);
            p += ((h & idx_mask) == b), s += ((h & idx_mask) != b);
        }

    return {p, s};
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool resize_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::lock(const std::size_t b)
{
    if constexpr(mt_ || resize_)
        if(mt_ || resize_worker_c > 1)
            while(true)
            {
                auto& w = (resize_ ? md_new(b) : md(b)).min_coord[0];
                const auto exp_unlocked = w & ~Metadata::lock_mask;
                if(__sync_bool_compare_and_swap(&w, exp_unlocked, exp_unlocked | Metadata::lock_mask))
                    break;
            }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool resize_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::lock_ordered(const std::size_t p, const std::size_t q)
{
    if constexpr(mt_ || resize_)
        KACHE_HASH_LIKELY(p != q) ?
            lock<resize_>(std::min(p, q)), lock<resize_>(std::max(p, q)) :
            lock<resize_>(p);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool resize_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::unlock(const std::size_t b)
{
    if constexpr(mt_ || resize_)
        if(mt_ || resize_worker_c > 1)
        {
            if constexpr(!resize_)
            {
                assert(md(b).min_coord[0] & Metadata::lock_mask);
                md(b).min_coord[0] &= ~Metadata::lock_mask;
            }
            else
            {
                assert(md_new(b).min_coord[0] & Metadata::lock_mask);
                md_new(b).min_coord[0] &= ~Metadata::lock_mask;
            }
        }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool resize_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::unlock_ordered(const std::size_t p, const std::size_t q)
{
    if constexpr(mt_ || resize_)
        KACHE_HASH_LIKELY(p != q) ?
                unlock<resize_>(std::max(p, q)), unlock<resize_>(std::min(p, q)) :
                unlock<resize_>(p);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::lock_sorted(std::size_t x, std::size_t y, std::size_t z)
{
    if constexpr(mt_)
    {
        if(x > y)   std::swap(x, y);
        if(y > z)   std::swap(y, z);
        if(x > y)   std::swap(x, y);

        lock(x);
        if(y != x)  lock(y);
        if(z != y)  lock(z);
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::unlock_sorted(const std::size_t x, const std::size_t y, const std::size_t z)
{
    if constexpr(mt_)
    {
        unlock(x);
        if(y != x)  unlock(y);
        if(z != x && z != y)    unlock(z);
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline std::size_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::size() const
{
    std::size_t sz = approx_sz;
    std::for_each(checkpoint.cbegin(), checkpoint.cend(), [&](auto& c){ sz += c.unwrap(); });
    return sz;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::find_in_overflow(const Kmer<k> key) const -> find_ret_t
{
    if constexpr(is_set_)
        return ov->find(key);
    else
    {
        const auto p = ov->find(key);
        return !p ? null_val : val_t(p->val);
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::place_at(const flat_t& e, uint8_t c, uint8_t m, std::size_t i, std::size_t j)
{
    slot(i, j) = e, md(i).cs[j] = c, md(i).min_coord[j] = m | (md(i).min_coord[j] & Metadata::lock_mask);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::place_at_resize(const flat_t& e, const uint8_t c, const uint8_t m, const std::size_t i, const std::size_t j)
{
    slot_new(i, j) = e, md_new(i).cs[j] = c, md_new(i).min_coord[j] = m | (md_new(i).min_coord[j] & Metadata::lock_mask);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::try_insert_at_resize(const flat_t& e, const uint8_t c, const uint8_t m, const std::size_t i)
{
    lock<true>(i);

    const auto c_vec = load_cs(md_new(i));
    const auto zero = zero_cs();
    const auto cmp = eq_mask(c_vec, zero);

    const auto j = std::countr_zero(cmp);
    if(j == B)
    {
        unlock<true>(i);
        return false;
    }

    place_at_resize(e, c, m, i, j);
    unlock<true>(i);
    return true;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::insert_at_resize(const flat_t& x, const uint8_t c, const uint8_t m)
{
    const auto idx_mask = cap_ - 1;
    const auto h = minimizer_hash(key(x), m);

    const auto b = h & idx_mask;
    if(try_insert_at_resize(x, c, m, b))
        return;

    const auto [s_1, s_2] = secondaries(h);
    const auto p = std::min(s_1, s_2), q = std::max(s_1, s_2);

    lock_ordered<true>(p, q);

    const auto zero = zero_cs();

    const auto c_vec_p = load_cs(md_new(p));
    const auto c_vec_q = load_cs(md_new(q));

    const auto cmp_p = eq_mask(c_vec_p, zero);
    const auto cmp_q = eq_mask(c_vec_q, zero);

    const auto sz_p = static_cast<std::size_t>(std::countr_zero(cmp_p));
    const auto sz_q = static_cast<std::size_t>(std::countr_zero(cmp_q));
    assert(std::countr_zero(cmp_p) + std::countl_one(cmp_p) == B), assert(std::countr_zero(cmp_q) + std::countl_one(cmp_q) == B);

    if(sz_p < sz_q)
    {
        place_at_resize(x, c, m, p, sz_p);
        unlock_ordered<true>(p, q);
    }
    else if(sz_q < B)
    {
        place_at_resize(x, c, m, q, sz_q);
        unlock_ordered<true>(p, q);
    }
    else
    {
        unlock_ordered<true>(p, q);
        if constexpr(is_set_)
            ov_new->insert(x, Overflow_Set_Val(c, m));
        else
            ov_new->insert(x.first, Overflow_Map_Val(c, m, x.second));
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::add_key(const uint16_t id)
{
    auto& cp = checkpoint[id].unwrap();
    if(KACHE_HASH_LIKELY(++cp < checkpoint_step.load(std::memory_order_relaxed)))
        return false;

    const auto sz = (approx_sz += cp);
    cp = 0;
    return sz >= resize_th.load(std::memory_order_relaxed);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::try_resize()
{
    if constexpr(mt_)   table_lock.lock();

    if(approx_sz >= resize_th)
        resize();

    if constexpr(mt_)   table_lock.unlock();
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline uint64_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::minimizer_hash(const Kmer<k> key, const uint8_t m)
{
    const uint8_t min_idx = m & coord_idx_mask;
    const bool in_fwd = m & coord_fwd_mask;
    const auto min = in_fwd ?
                        key.template lmer_at<l>(min_idx) : Kmer_Utility::reverse_complement<l>(key.template lmer_at<l>(k - min_idx - l));
    return wang_u64(min);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::resize()
{
    const auto old_cap = cap_;
    cap_ *= 2;
    T_new = aligned_allocate<flat_t>(cap_ * B, L1_CACHE_LINE_SIZE);
    M_new = aligned_allocate<Metadata>(cap_, L1_CACHE_LINE_SIZE);
    ov_new = new ov_t(capacity() * of_default);
    std::memset(M_new, 0, cap_ * sizeof(Metadata)); // Checksum `0` denotes absence of a key.

    std::vector<std::thread> W;
    W.reserve(resize_worker_c);
    for(uint64_t i = 0; i < resize_worker_c; ++i)
        W.emplace_back([&](const auto w_id)
        {
            const auto stride = resize_worker_c;
            for(std::size_t b = w_id; b < old_cap; b += stride)
                for(std::size_t j = 0; j < B; ++j)
                    if(md(b).cs[j])
                        insert_at_resize(slot(b, j), md(b).cs[j], md(b).min_coord[j] & ~Metadata::lock_mask);

            auto it = ov->iterator(resize_worker_c, w_id);
            typedef std::conditional_t<!is_map_, Overflow_Set_Val, Overflow_Map_Val> ov_val_t;
            Kmer<k> key;
            ov_val_t val;
            while(it.next(key, val))
                if constexpr(is_set_)
                    insert_at_resize(key, val.cs, val.min_coord);
                else
                    insert_at_resize(flat_t(key, val.val), val.cs, val.min_coord);
        }, i);

    std::for_each(W.begin(), W.end(), [](auto& w){ w.join(); });


    std::swap(M, M_new), std::swap(T, T_new), std::swap(ov, ov_new);
    deallocate(M_new), deallocate(T_new), delete(ov_new);
    M_new = nullptr, T_new = nullptr, ov_new = nullptr;

    resize_th = capacity() * lf;
    checkpoint_step = step_for(resize_th);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename Op_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::stream(const char* const* const seqs, const std::size_t* const lens, const std::size_t c, const Op_& op)
{
    for(std::size_t q = 0; q < c; ++q)
    {
        const char* const seq = seqs[q];
        const auto len = lens[q];
        for(std::size_t s = 0; s + k <= len; )  // `s` is the start of the current stretch of ACGT bases.
        {
            const auto e = s + DNA_Utility::acgt_prefix(seq + s, len - s);   // End of the stretch.

            if(e - s >= k)
            {
                Kmer_Batch<k, l> kb;
                typename Kmer_Batch<k, l>::Chunk ch;
                kb.reset(seq + s, e - s - (k - 1));
                for(std::size_t i = s; kb.next(ch); i += ch.n)
                    for(std::size_t j = 0; j < ch.n; ++j)
                        op(Kmer_Desc{Kmer<k>(ch.key[j]), ch.h[j], ch.c[j], ch.m[j]}, q, i + j);
            }

            s = e + 1;
        }
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline uint32_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::probe(const Kmer<k> key, const uint8_t c, const std::size_t i, const cs_vec c_vec, uint32_t& e) const
{
    assert(i < cap_);

    auto cmp = eq_mask(c_vec, splat_cs(c));
    e = std::countr_zero(eq_mask(c_vec, zero_cs()));

    while(cmp)
    {
        const auto j = std::countr_zero(cmp);
        cmp &= (cmp - 1);
        if(key_equals(slot(i, j), key))
            return j;
    }

    return B;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
std::pair<std::size_t, std::size_t> Streaming_Kmer_Hash_Table<k, mt_, T_, l>::secondaries(const uint64_t h) const
{
    const auto idx_mask = cap_ - 1;
    return {XXH3_64bits_withSeed(&h, sizeof(h), 0) & idx_mask, XXH3_64bits_withSeed(&h, sizeof(h), 1llu << 63) & idx_mask};
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename F, typename V>
bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::ov_update(const Kmer<k> key, const uint8_t c, const uint8_t m, const F& f, const V& val)
{
    if constexpr(is_set_)
        return (void)f, (void)val, !ov->insert(key, Overflow_Set_Val(c, m));
    else
        return !ov->upsert(key, Overflow_Map_Val(c, m, val), [&](const Overflow_Map_Val& v){ return Overflow_Map_Val(v.cs, v.min_coord, f(v.val)); });
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename F, typename V>
__attribute__((flatten)) inline std::size_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::update(const char* const* const seqs, const std::size_t* const lens, const std::size_t c, const F& f, const V& val, const Token& token)
{
    // The current minimizer run, its primary bucket `b` and secondary buckets `p` and `q`, and the buckets this
    // thread holds locked for it: none (0), `b` (1), or all three (3). Once `b` is found full, the run's remaining
    // k-mers go through all three.
    uint64_t run_h = 0;
    bool in_run = false;
    bool spilled = false;
    std::size_t b = 0, p = 0, q = 0;
    int held = 0;
    cs_vec cs_b;    // The checksums of `b` while this thread holds it alone, kept in registers. Otherwise a single
                    // checksum's store followed by the next probe's load of all of them may block store forwarding.

    const auto release = [&]()
    {
        if(held == 1)
            unlock(b);
        else if(held == 3)
            unlock_sorted(b, p, q);

        held = 0;
    };

    std::size_t ins_c = 0;  // Number of k-mers inserted.

    // Resizes the table, iff it is due.
    const auto resize_now = [&]()
    {
        release();
        if constexpr(mt_)   table_lock.unlock_shared(token.id);
        try_resize();
        if constexpr(mt_)   table_lock.lock_shared(token.id);
        in_run = false;
    };

    // Accounts for a k-mer inserted, resizing the table if it is due.
    const auto inserted = [&]()
    {
        ins_c++;
        if(KACHE_HASH_UNLIKELY(add_key(token.id)))
            resize_now();
    };

    // Places the k-mer `d` into the `i`'th bucket at slot `j`.
    const auto place = [&](const Kmer_Desc& d, const std::size_t i, const std::size_t j)
    {
        if constexpr(is_set_)
            place_at(d.key, d.c, d.m, i, j);
        else
            place_at(flat_t(d.key, val), d.c, d.m, i, j);
    };

    // Applies the update to the value of the present k-mer at the `i`'th bucket's slot `j`.
    const auto found = [&](const std::size_t i, const std::size_t j)
    {
        if constexpr(is_map_)
            slot(i, j).second = f(slot(i, j).second);
        else
            (void)i, (void)j;
    };

    // Applies the overflow table's update with the k-mer `d`, holding no lock.
    const auto apply_overflow = [&](const Kmer_Desc& d)
    {
        if(ov_update(d.key, d.c, d.m, f, val))
            inserted();
    };

    // Applies the update with the k-mer `d`.
    const auto apply = [&](const Kmer_Desc& d)
    {
        if(!in_run || d.h != run_h)
        {
            release();
            in_run = true, run_h = d.h, spilled = false;
            b = d.h & (cap_ - 1);
        }

        if(!spilled)
        {
            if(held == 0)
                lock(b), held = 1, cs_b = load_cs(md(b));

            uint32_t e;
            if(const auto j = probe(d.key, d.c, b, cs_b, e); j < B)
                return found(b, j);

            if(e < B)
            {
                place(d, b, e);
                cs_b = set_cs(cs_b, e, d.c);
                return inserted();
            }

            unlock(b), held = 0;
            std::tie(p, q) = secondaries(d.h);
            spilled = true;
        }

        if(held == 0)
            lock_sorted(b, p, q), held = 3;

        uint32_t e_b, e_p, e_q;
        if(const auto j = probe(d.key, d.c, b, e_b); j < B)
            return found(b, j);
        if(const auto j = probe(d.key, d.c, p, e_p); j < B)
            return found(p, j);
        if(const auto j = probe(d.key, d.c, q, e_q); j < B)
            return found(q, j);

        if(e_p < e_q)
            return place(d, p, e_p), inserted();
        if(e_q < B)
            return place(d, q, e_q), inserted();

        release();
        apply_overflow(d);
    };

    if constexpr(mt_)   table_lock.lock_shared(token.id);

    stream(seqs, lens, c, [&](const Kmer_Desc& d, std::size_t, std::size_t)
    {
        if constexpr(mt_)
            if(KACHE_HASH_UNLIKELY(table_lock.writer_pending()))  // Yield to a resize.
            {
                release();
                table_lock.unlock_shared(token.id);
                table_lock.lock_shared(token.id);
                in_run = false;
            }

        apply(d);
    });

    release();
    if constexpr(mt_)   table_lock.unlock_shared(token.id);

    return ins_c;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename G>
__attribute__((flatten)) inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::find(const char* const* const seqs, const std::size_t* const lens, const std::size_t c, const G& g) const
{
    // The current minimizer run, its primary bucket `b`, and its secondary buckets `p` and `q` once computed.
    uint64_t run_h = 0;
    bool in_run = false;
    bool sec = false;
    std::size_t b = 0, p = 0, q = 0;

    // Returns what `find` returns for the present k-mer at the `i`'th bucket's slot `j`.
    const auto present = [&](const std::size_t i, const std::size_t j) -> find_ret_t
    {
        (void)i, (void)j;
        if constexpr(is_map_)   return slot(i, j).second;
        else                    return true;
    };

    stream(seqs, lens, c, [&](const Kmer_Desc& d, const std::size_t seq_i, const std::size_t i)
    {
        if(!in_run || d.h != run_h)
            in_run = true, run_h = d.h, sec = false, b = d.h & (cap_ - 1);

        uint32_t e;
        if(const auto j = probe(d.key, d.c, b, e); j < B)
            return g(seq_i, i, present(b, j));
        if(e < B)
            return g(seq_i, i, find_ret_t(null_val));

        if(!sec)
            std::tie(p, q) = secondaries(d.h), sec = true;

        uint32_t e_p, e_q;
        if(const auto j = probe(d.key, d.c, p, e_p); j < B)
            return g(seq_i, i, present(p, j));
        if(const auto j = probe(d.key, d.c, q, e_q); j < B)
            return g(seq_i, i, present(q, j));

        g(seq_i, i, (e_p == B && e_q == B) ? find_in_overflow(d.key) : find_ret_t(null_val));
    });
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline std::size_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::count(const char* const* const seqs, const std::size_t* const lens, const std::size_t c) const
{
    std::size_t hit_c = 0;
    find(seqs, lens, c, [&](std::size_t, std::size_t, const find_ret_t& r){ hit_c += static_cast<bool>(r); });
    return hit_c;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::key_distribution() const -> std::tuple<double, double, double>
{
    std::size_t p = 0, s = 0;
    for(size_t b = 0; b < cap_; ++b)
    {
        const auto d = distribution(b);
        p += d.first, s += d.second;
    }

    const auto o = overflow_size();
    const auto sz = p + s + o;
    return {p * 100.0 / sz, s * 100.0 / sz, o * 100.0 / sz};
}


}



#endif
