
#ifndef STREAMING_KMER_HASHSET_HPP
#define STREAMING_KMER_HASHSET_HPP



#include "DNA_Utility.hpp"
#include "Kmer.hpp"
#include "Kmer_Batch.hpp"
#include "Kmer_Lanes.hpp"
#include "Concurrent_Hash_Table.hpp"
#include "RW_Lock.hpp"
#include "utility.hpp"

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <bit>
#include <atomic>
#include <tuple>
#include <optional>
#include <utility>
#include <type_traits>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <thread>
#include <algorithm>
#include <barrier>
#include <cassert>
#include <cmath>

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
template <uint16_t k, bool mt_, typename T_ = void, uint16_t l = 22>
class Streaming_Kmer_Hash_Table
{
    static_assert(k <= 32, "The k-mers need to fit in a machine word.");

public:

    class Token;

    static constexpr uint64_t max_user = 128;   // Maximum number of users allowed to use the table: a thread per core of a 2×64-core node, for now.
    static constexpr double lf_default_ = 0.8;  // Default load factor: 80%.

private:

    static constexpr double of_default = 1.0 / 32;  // Initial overflow factor: the overflow table starts with room for
                                                    // this fraction of the main table's capacity, and grows on its own.
    static constexpr double ov_lf = 0.75;   // Load factor at which the overflow table grows.
    static constexpr std::size_t ov_min_cap = 1024; // Fewest keys the overflow table has room for.
    static constexpr std::size_t B = 32;    // Bucket size: a bucket's checksums fill a 256-bit vector.

    typedef uint32_t mask_t;    // A bit per slot of a bucket.

    static constexpr bool is_set_ = std::is_void_v<T_>; // Whether a hash-set.
    static constexpr bool is_map_ = !is_set_;   // Whether a hash-map.
    typedef std::conditional_t<is_map_, std::pair<Kmer<k>, T_>, Kmer<k>> flat_t;    // An element of the table: a key, with its value in a map.
    typedef std::optional<T_> val_t;    // Type of values returned by hash-map operations.
    typedef std::conditional_t<is_map_, val_t, bool> find_ret_t;    // Return type for find methods.
    static constexpr auto null_val = [](){ if constexpr(is_map_) return std::nullopt;   else return false; }(); // Marker for absent keys (also values).

    static constexpr auto& key(const flat_t& e) { if constexpr(is_map_) return e.first; else return e; }    // Returns the key from a table-element `e`.

    const double lf;    // Load factor.

    std::size_t cap_;   // Capacity of the table in terms of buckets.
    static constexpr std::size_t max_cap = (std::size_t(1) << 32) - 1;  // Most buckets of a table: range reduction
                                                                        // takes 32-bit products.
    std::atomic_uint64_t resize_th; // At what size of the table to trigger the next resize.

    // Metadata of the k-mers in a bucket: their checksums, and 256 flag bits.
    // Flag 0 locks the bucket, and the others are its spill filter: they mark
    // the minimizers whose keys left it, their primary bucket.
    struct Metadata
    {
        uint8_t cs[B];  // Checksums of the k-mers.
        uint8_t flags[32];  // Flag `i` is bit `i % 8` of `flags[i / 8]`.

        static constexpr uint8_t lock_mask = 0b0000'0001;   // Bit of `flags[0]` used to lock the bucket.
    };

    struct Overflow_Set_Val {}; // A set's overflow table maps its keys to nothing.
    typedef std::conditional_t<is_set_, Overflow_Set_Val, T_> ov_val_t; // Value type of the overflow table; it keeps
                                                                        // no checksums, which resizes recompute.
    typedef Concurrent_Hash_Table<Kmer<k>, ov_val_t, Hash<Kmer<k>>> ov_t;   // The overflow table type.


    typedef std::conditional_t<is_map_, T_, uint8_t> mapped_t;  // The mapped type; a placeholder in a set.
    static_assert(is_set_ || std::is_trivially_copyable_v<mapped_t>, "Values may be copied bytewise into and out of the slots.");

    // Storage of a table's buckets: each bucket holds the metadata of its
    // keys, then their slots, so that an operation finds both on one page. A
    // slot holds a key and, in a map, its value. Where the elements would pad,
    // as with values narrower than keys, the slots pack them bytewise, so that
    // a value takes only its own bytes. The buckets sit on huge pages where
    // supported.
    class Store
    {
        static constexpr std::size_t line = L1_CACHE_LINE_SIZE;
        static constexpr std::size_t key_sz = sizeof(Kmer<k>);  // Bytes of a key.
        static constexpr std::size_t val_sz = (is_map_ ? sizeof(mapped_t) : 0);    // Bytes of a value.
        static constexpr std::size_t slot_sz = key_sz + val_sz; // Bytes of a slot.
        static constexpr bool packed = (sizeof(flat_t) > slot_sz);  // Whether the slots are packed bytewise.

        // A bucket: the metadata of its keys, then their slots.
        struct alignas(L1_CACHE_LINE_SIZE) Bucket
        {
            Metadata md;    // Metadata of the keys.
            std::conditional_t<packed, unsigned char[B * slot_sz], flat_t[B]> slot; // The keys, with their values in a map.
        };

        static constexpr std::size_t line_slots = std::max<std::size_t>(1, line / sizeof(flat_t)); // Number of unpacked slots in a cache line.

        // The packed slots whose insertions prefetch the line after their own,
        // as a mask: each line's first slot to start in it, if a line follows
        // within the slots.
        static constexpr mask_t line_starts = []()
        {
            mask_t m = 0;
            for(std::size_t e = 0; e < B; ++e)
                if((e * slot_sz) % line < slot_sz && ((e * slot_sz) / line + 1) * line < B * slot_sz)
                    m |= mask_t(1) << e;
            return m;
        }();

        // The packed slots before which the first insertion of a run
        // prefetches the line after its own too.
        static constexpr std::size_t run_lim = []()
        {
            std::size_t e = 0;
            while(e < B && e * slot_sz + line < B * slot_sz)
                e++;
            return e;
        }();

        Bucket* T = nullptr;    // The buckets.

        // Returns the `j`'th packed slot of the `i`'th bucket.
        unsigned char* at(std::size_t i, std::size_t j) const;

    public:

        static constexpr std::size_t bucket_bytes = sizeof(Bucket); // Bytes of a bucket.

        // Returns the metadata of the `i`'th bucket.
        Metadata& md(const std::size_t i) const { return T[i].md; }

        // Returns the key at the `j`'th slot of the `i`'th bucket.
        Kmer<k> key(std::size_t i, std::size_t j) const;

        // Returns the value at the `j`'th slot of the `i`'th bucket, in a map.
        mapped_t val(std::size_t i, std::size_t j) const;

        // Sets the value at the `j`'th slot of the `i`'th bucket to `v`, in a
        // map.
        void set_val(std::size_t i, std::size_t j, const mapped_t& v) const;

        // Sets the value `v` at the `j`'th slot of the `i`'th bucket to
        // `f(v)`, in a map, as `update_in_place` does.
        template <typename F> void update_val(std::size_t i, std::size_t j, const F& f) const;

        // Returns the element at the `j`'th slot of the `i`'th bucket.
        flat_t get(std::size_t i, std::size_t j) const;

        // Puts the element `e` at the `j`'th slot of the `i`'th bucket.
        void put(std::size_t i, std::size_t j, const flat_t& e) const;

        // Prefetches the `j`'th slot of the `i`'th bucket, with write intent
        // iff `write_`.
        template <bool write_> void prefetch_slot(std::size_t i, std::size_t j) const;

        // Prefetches what insertions from the `e`'th slot of the `i`'th bucket
        // touch first: the slot's line, and the one after.
        void prefetch_run(std::size_t i, std::size_t e) const;

        // Prefetches the line that a run of insertions reaches next once it
        // has placed a key at the `e`'th slot of the `i`'th bucket, iff the
        // slot is the first one to start in its line.
        void prefetch_frontier(std::size_t i, std::size_t e) const;

        // Allocates the storage of `n` buckets, zeroed, and maps it in
        // parallel with `t_c` threads.
        void allocate(std::size_t n, std::size_t t_c);

        // Grows the storage from `n` to `m` buckets, keeping the first `n`
        // buckets; the new ones come zeroed. Returns whether it succeeded.
        bool grow(std::size_t n, std::size_t m);

        // Deallocates the storage of `n` buckets.
        void deallocate(std::size_t n);
    };

    Store T;    // The table's `cap_` buckets.
    ov_t* ov;   // The overflow table.

    Store T_new;    // New buckets after a resize.
    ov_t* ov_new;    // New overflow table after a resize.

    // Returns the metadata of the `i`'th bucket.
    Metadata& md(const std::size_t i) const { return T.md(i); }

    // Returns the metadata of the `i`'th new bucket in a resize.
    Metadata& md_new(const std::size_t i) const { return T_new.md(i); }

    std::atomic_uint16_t registered_thread_count;   // Number of threads registered to use the table.

    std::atomic_uint64_t approx_sz; // Size of the hash table, short of the user-local counts not yet flushed into it.
    std::vector<Padded<uint64_t>> checkpoint;   // `checkpoint[i]` is the number of keys the `i`'th user has added and
                                                // not yet flushed into `approx_sz`, as of its last update's end: an
                                                // update counts in registers.
    static constexpr uint64_t resize_checkpoint = 16384;    // Most keys a user adds before flushing its count and checking for a resize.
    std::atomic_uint64_t checkpoint_step;   // Number of keys a user adds before flushing its count and checking for a resize: unflushed keys stay under an eighth of the resize threshold.
    RW_Lock<max_user> table_lock; // Readers-writers lock to resize the table exclusively.

    const uint64_t resize_worker_c; // Number of worker threads to spawn for a concurrent resize.

    std::atomic_uint64_t ov_sz; // Size of the overflow table, short of the user-local counts not yet flushed into it.
    std::vector<Padded<uint64_t>> ov_checkpoint;    // `ov_checkpoint[i]` is the number of keys the `i`'th user has added to the overflow table and not yet flushed into `ov_sz`.
    std::atomic_uint64_t ov_th; // Size of the overflow table at which it grows.
    std::atomic_uint64_t ov_step;   // Number of keys a user adds to the overflow table before flushing its count and checking for growth.


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

    // Returns the number of buckets of a table for `max_sz` keys at load
    // factor `lf`: the fewest whose resize threshold exceeds `max_sz`, rounded
    // up to a power of 2 iff `pow2`.
    static std::size_t cap_for(std::size_t max_sz, double lf, bool pow2);

    // Returns the primary bucket of the minimizer hash `h` in a table of `cap`
    // buckets: the high half of `cap` times the hash's low 32 bits. Doubling
    // `cap` maps the bucket `b` to `2b` or `2b + 1`.
    // NB: Not the hash's top bits, which window minima bias.
    static std::size_t bucket_of(uint64_t h, std::size_t cap);

    // Returns the primary bucket of the minimizer hash `h`.
    std::size_t bucket_of(uint64_t h) const;

    // Puts the primary buckets of the minimizer hashes `h[0..n)` in a table of
    // `cap` buckets at `b[0..n)`.
    static void buckets_of(const uint64_t* h, std::size_t n, std::size_t cap, std::size_t* b);

    // Returns the primary buckets of the minimizer hashes in the lanes `h`, in
    // a table of `cap` buckets.
    static typename Kmer_Lanes<k, l>::vec buckets_of(typename Kmer_Lanes<k, l>::vec h, std::size_t cap);

    // Returns the bucket that the uniform 64-bit hash `x` picks in a table of
    // `cap` buckets: the high word of their product.
    static std::size_t pick_bucket(uint64_t x, std::size_t cap);

    // Returns the index in `[1, 256)` of the flag that marks the minimizer hash
    // `h` in a spill filter.
    // NB: From bits 32–39 of `h`: the bucket takes its low 32 bits, and window
    // minima bias its top ones.
    static std::size_t filter_pos(uint64_t h);

    // Returns whether the spill filter of the bucket with metadata `md` marks
    // the minimizer hash `h`.
    static bool may_spill(const Metadata& md, uint64_t h);

    // Marks the minimizer hash `h` in the spill filter of the bucket with
    // metadata `md`.
    static void mark_spill(Metadata& md, uint64_t h);

    // Clears the spill filter of the bucket with metadata `md`.
    static void clear_spills(Metadata& md);

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

    // Sets the value `v` to `f(v)`, in a map: atomically in a concurrent
    // table, where `f` may be applied more than once.
    template <typename F> static void update_in_place(mapped_t& v, const F& f);

    // Sets the associated value `v` of the k-mer `key` to `f(v)` iff it is
    // present in the overflow table, as `update_in_place` does. Returns whether
    // it is present.
    // NB: Out of line, as the streaming updates rarely reach the overflow
    // table.
    template <typename F> __attribute__((noinline)) bool update_in_overflow(Kmer<k> key, const F& f);

    // Directly places the element `e` with checksum `c` into the `i`'th bucket
    // at index `j`.
    void place_at(const flat_t& e, uint8_t c, std::size_t i, std::size_t j);

    // Directly places the element `e` with checksum `c` into the `i`'th bucket
    // at index `j` of the new resized table.
    void place_at_resize(const flat_t& e, uint8_t c, std::size_t i, std::size_t j);

    // Tries to insert the element `x` with checksum `c` and minimizer hash `h`
    // into its primary bucket `i` of the newly resized table. Returns whether
    // there was empty space; if not, marks `h` in the bucket's spill filter.
    bool try_insert_at_resize(const flat_t& x, uint8_t c, uint64_t h, std::size_t i);

    // Inserts the element `x` with checksum `c` and minimizer hash `h` into
    // the resized table, whose secondary buckets are `p` and `q` iff `sec`, as
    // computed ahead. Returns whether it goes to the overflow table.
    // NB: In a resize, buckets only gain keys, so one found full without its
    // lock stays so: overflow-bound keys take no bucket locks.
    bool insert_at_resize(const flat_t& x, uint8_t c, uint64_t h, bool sec, std::size_t p, std::size_t q);

    // Puts the minimizer hashes of the `n` keys at the `b`'th bucket at
    // `h[0..n)`, recomputed from the keys.
    void bucket_hashes(std::size_t b, std::size_t n, uint64_t* h) const;

    // Returns the secondary buckets of the minimizer hash `h`.
    __attribute__((noinline)) std::pair<std::size_t, std::size_t> secondaries(uint64_t h) const;

    // Returns the number of keys a user adds before flushing its count and
    // checking for a resize, for the resize threshold `th`.
    static uint64_t step_for(uint64_t th) { return std::clamp<uint64_t>(th / (8 * max_user), 1, resize_checkpoint); }

    // Counts a key newly added to the overflow table by the `id`'th user.
    // Returns whether the overflow table is due to grow.
    bool add_ov_key(uint16_t id);

    // Sets the overflow table's growth threshold and checkpoint step, for its
    // current capacity.
    void set_ov_th();

    // Returns the number of keys to size a resize's new overflow table for: as
    // many as the current one, as most of its keys stay overflowed; and at
    // least the overflow factor's share of the new capacity.
    std::size_t new_ov_size() const;

    // Resets the overflow table's size to `sz`, after a resize: none of the
    // users' counts are unflushed then.
    void reset_ov_size(uint64_t sz);

    // Doubles the overflow table iff it is due, or iff it is still `ov_seen`
    // and `full`, i.e. found full by an insertion.
    __attribute__((noinline, cold)) void try_grow_ov(const ov_t* ov_seen, bool full);

    // Tries to resize the hash table to double its capacity: does so iff it is
    // due.
    __attribute__((noinline, cold)) void try_resize();

    // Resizes the hash table to double its capacity in place: the arrays grow,
    // and a key at its primary bucket `b` moves to `2b` or `2b + 1`; else it is
    // set aside and re-inserted.
    void resize();

    // Descriptor of a k-mer for the streaming operations: what the table
    // operations need of it.
    struct Kmer_Desc
    {
        Kmer<k> key;    // Canonical form of the k-mer.
        uint64_t h; // Hash of its minimizer.
        uint8_t c;  // Its checksum.
        std::size_t b;  // Its primary bucket at the stream's capacity, iff the stream computes them.
    };

    // Returns the primary bucket of the k-mer `d`, as described.
    static std::size_t bucket_of(const Kmer_Desc& d);

    // Whether the streaming operation `Op_` takes the k-mers' orientations, as
    // `op(d, q, i, fwd)`.
    template <typename Op_>
    static constexpr bool oriented_op = std::is_invocable_v<const Op_&, const Kmer_Desc&, std::size_t, std::size_t, bool>;

    static constexpr std::size_t slot_dist = 16;    // How many k-mers ahead of their table operations the streaming
                                                    // operations prefetch the slots of their keys in their primary
                                                    // buckets, or else what they touch beyond. In the lanes, that
                                                    // prefetches all of a chunk's slots a chunk ahead.

    // Applies `op(d, q, i)` to the descriptor `d` of each k-mer of the
    // sequences `seqs[0..c)` of lengths `lens[0..c)`, the k-mer being at index
    // `i` of the `q`'th sequence; or `op(d, q, i, fwd)` iff `op` takes it, with
    // `fwd` as whether the k-mer, as it occurs, is its canonical form, i.e.
    // its key. Non-ACGT characters split the sequences. A
    // stretch's k-mers come in order, but stretches may interleave. Applies
    // `pre_meta(d)` to each descriptor some k-mers ahead of `op`, and
    // `pre_slot(make)` `slot_dist` k-mers ahead, where `make()` returns the
    // descriptor. Skipping a prefetch then skips making its descriptor. Iff
    // `buckets_`, the descriptors carry their primary buckets in a table of
    // `cap` buckets. With `lanes_`, the k-mers are described in SIMD lanes, and
    // iff `cont_`, lane runs continue across chunks; else the stretches of ACGT
    // bases go through one pipeline of the batch describer. If `pre_slot` has
    // `active()`, the lanes skip a chunk's slot prefetches while it is `false`.
    // NB: `pre_meta` takes every descriptor, not just each minimizer run's
    // first: a prefetch per k-mer costs less than branching on run starts.
    template <bool lanes_, bool cont_, bool buckets_, typename Pre_Meta_, typename Pre_Slot_, typename Op_>
    static void stream(const char* const* seqs, const std::size_t* lens, std::size_t c, std::size_t cap, const Pre_Meta_& pre_meta, const Pre_Slot_& pre_slot, const Op_& op);

    // Applies `stream`'s callbacks to the k-mers of the sequences `seqs[0..c)`
    // of lengths `lens[0..c)` through the SIMD lanes: long stretches split over
    // the lanes on their own, short ones gathered a lane each.
    template <bool cont_, bool buckets_, typename Pre_Meta_, typename Pre_Slot_, typename Op_>
    static void stream_lanes_all(const char* const* seqs, const std::size_t* lens, std::size_t c, std::size_t cap, const Pre_Meta_& pre_meta, const Pre_Slot_& pre_slot, const Op_& op);

    // Applies `stream`'s callbacks to the `n` k-mers of the ACGT stretch `t`,
    // whose first k-mer has index `base` in the `q`'th sequence, describing
    // them in chunks: the metadata of a chunk's buckets is prefetched one to
    // two chunks ahead.
    template <bool buckets_, typename Pre_Meta_, typename Pre_Slot_, typename Op_>
    static void stream_batched(const char* t, std::size_t n, std::size_t q, std::size_t base, std::size_t cap, const Pre_Meta_& pre_meta, const Pre_Slot_& pre_slot, const Op_& op);

    static constexpr std::size_t lane_c = Kmer_Lanes<k, l>::L;  // Number of SIMD lanes describing stretches at once.
    static constexpr std::size_t lane_min_kmers = 4096; // Fewest k-mers in a stretch for it to be split over the lanes on its own.
    static constexpr std::size_t lane_depth = 2;    // Number of lane chunks in flight: the one operated on, and the one
                                                    // ahead, whose metadata is prefetched once described.
    static constexpr std::size_t lane_group_c = 16; // Most groups of short stretches, such as reads, that go through one
                                                    // pipeline of lane chunks: 1 restarts the pipeline per group.

    // The streaming operations describe k-mers in SIMD lanes only for targets
    // with AVX-512, whose cores run 256-bit vectors at full width. Cores that
    // split them into 128-bit halves, such as Sierra Forest's E-cores, run the
    // batch describer faster. Tests may force either with `KACHE_HASH_LANES`.
#if defined(KACHE_HASH_LANES)
    static constexpr bool use_lanes = (KACHE_HASH_LANES != 0);
#elif defined(__AVX512F__)
    static constexpr bool use_lanes = true;
#else
    static constexpr bool use_lanes = false;
#endif

    // A stretch of ACGT bases: at `t`, of `len` bases, and its first k-mer at
    // index `base` of the `q`'th sequence.
    struct Stretch
    {
        const char* t;
        std::size_t len;
        std::size_t q;
        std::size_t base;
    };

    // Applies `stream`'s callbacks to the k-mers of the stretches `S[0..c)` in
    // order, describing them in chunks as `stream_batched` does, but in one
    // pipeline: the next stretch is described, and its metadata prefetched,
    // while operating on the current one.
    template <bool buckets_, typename Pre_Meta_, typename Pre_Slot_, typename Op_>
    static void stream_stretches(const Stretch* S, std::size_t c, std::size_t cap, const Pre_Meta_& pre_meta, const Pre_Slot_& pre_slot, const Op_& op);

    // Applies `stream`'s callbacks to the k-mers of `g_c` groups of
    // stretches: the `g`'th group's are `S[g * lane_c + j]` for `j` in
    // `[0, G_c[g])`, with `G_c[g] <= lane_c`. Each group goes through the SIMD
    // lanes at once, a stretch per lane, and the groups through one pipeline of
    // chunks. The metadata of a chunk's buckets is prefetched once it is
    // described, before the operations on the chunk ahead of it. A chunk's
    // k-mers go lane by lane, a lane's in order. All of a chunk's slots are
    // prefetched a chunk ahead. Iff `cont_`, lane runs continue across chunks.
    template <bool cont_, bool buckets_, typename Pre_Meta_, typename Pre_Slot_, typename Op_>
    static void stream_lanes(const Stretch* S, const std::size_t* G_c, std::size_t g_c, std::size_t cap, const Pre_Meta_& pre_meta, const Pre_Slot_& pre_slot, const Op_& op);

    // Describes the next chunk of the lanes `kl` into `out`, with the k-mers'
    // orientations iff `orient_`. Returns what `Kmer_Lanes::next` does.
    // NB: Out of line, so that the lanes' unrolled body stays apart from the
    // operations' loops. The batch describer runs faster inlined.
    template <bool orient_ = false>
    __attribute__((noinline)) static bool describe_lanes(Kmer_Lanes<k, l>& kl, typename Kmer_Lanes<k, l>::Chunk& out);

    // Returns the index of the slot holding the k-mer `key` with checksum `c`
    // in the `i`'th bucket, or `B` if it is absent there. Puts the index of the
    // bucket's first empty slot at `e`, `B` if it is full.
    uint32_t probe(Kmer<k> key, uint8_t c, std::size_t i, uint32_t& e) const { return probe(key, c, i, load_cs(md(i)), e); }

    // Returns the index of the slot holding the k-mer `key` with checksum `c`
    // in the `i`'th bucket, which has the checksums `c_vec`; or `B` if it is
    // absent there. Puts the index of the bucket's first empty slot at `e`, `B`
    // if it is full.
    uint32_t probe(Kmer<k> key, uint8_t c, std::size_t i, cs_vec c_vec, uint32_t& e) const;

    // Returns the index of the slot holding the k-mer `key` with checksum `c`
    // in the `i`'th bucket, or `B` if it is absent there; then sets `full` iff
    // the bucket is full, from its last checksum alone, as keys fill a bucket's
    // slots in order.
    uint32_t probe(Kmer<k> key, uint8_t c, std::size_t i, bool& full) const;

    // Prefetches the metadata of the bucket of the k-mer `d`, with write
    // intent iff `write_`, i.e. for operations that may insert. These compute
    // the bucket afresh, as their resizes change the capacity; others take
    // the one described.
    template <bool write_> void prefetch_meta(const Kmer_Desc& d) const;

    // Secondary buckets of the latest minimizer run seen by a prefetch stage.
    struct Sec_Cache
    {
        uint64_t h = 0; // Minimizer hash of the run.
        std::size_t cap = 0;    // Table capacity, in buckets, when they were computed; 0 if never.
        std::size_t p, q;   // The secondary buckets.
    };

    // Returns the secondary buckets of the minimizer hash `h`, through the
    // cache `sc`.
    std::pair<std::size_t, std::size_t> secondaries(uint64_t h, Sec_Cache& sc) const;

    // Prefetches the slot where an operation on the k-mer `d` likely ends in
    // its primary bucket: its key's, by checksum, or else the first empty one
    // iff it may insert, i.e. iff `insert_`; with write intent iff `write_`.
    // As in `prefetch_meta`, operations that may insert compute the bucket
    // afresh.
    // If the operation goes beyond the bucket, prefetches its secondary
    // buckets' metadata instead, computing them through the cache `sc`, and
    // its overflow-table slot. An operation that doesn't insert goes beyond
    // only if the spill filter marks its minimizer.
    template <bool write_, bool insert_ = write_> void prefetch_slot(const Kmer_Desc& d, Sec_Cache& sc) const;

    // Prefetches what an operation on the k-mer `key` with minimizer hash `h`
    // touches beyond its full primary bucket, as `prefetch_slot` does,
    // computing the secondary buckets through the cache `sc`.
    // NB: Out of line, as inlined it bloats the operations' loops by 12-18
    // instructions per k-mer. It takes the k-mer and the hash by value, as a
    // descriptor by reference would be kept in memory for the call.
    template <bool write_> __attribute__((noinline, cold)) void prefetch_spill(Kmer<k> key, uint64_t h, Sec_Cache& sc) const;

    // Lookups' slot prefetches, as `stream` takes them: on while the last
    // lookup hit, and with write intent iff `write_`, for updates in place.
    // The lanes ask once per chunk whether they are on.
    template <bool write_>
    struct Lookup_Pre_Slot
    {
        const Streaming_Kmer_Hash_Table& M; // The table.
        Sec_Cache& sc;  // The cache of secondary buckets of the prefetches.
        const bool& hit;    // Whether the last lookup hit.

        // Returns whether the prefetches are on.
        bool active() const;

        // Prefetches the slot of the k-mer that `make()` describes, as
        // `prefetch_slot` does, iff the prefetches are on.
        template <typename Make_> void operator()(const Make_& make) const;
    };

    // Inserts the k-mer `key` into the overflow table iff it is absent there,
    // with value `val` in a map; also in a map, updates the value `v` of a
    // present one to `f(v)`. Returns whether it was absent. Sets `full` iff it
    // is absent and finds no room; then nothing changes.
    // NB: Out of line, as the streaming updates rarely reach the overflow
    // table.
    template <typename F, typename V>
    __attribute__((noinline)) bool ov_update(Kmer<k> key, const F& f, const V& val, bool& full);

    // Updates the table with the k-mers of the sequences `seqs[0..c)` of
    // lengths `lens[0..c)`: inserts each absent k-mer with value `val`, and in
    // a map, updates the value `v` of each present one to `f(v)`. Iff
    // `orient_`, in a map, updates each k-mer's value `v` to `f(v, fwd)`
    // instead, from `val` for an absent one, with `fwd` as whether it occurs in
    // its canonical form. Returns the number of k-mers inserted. The invoking
    // user thread must provide its token `token`.
    // NB: The streaming operations force-inline all their callees but the cold
    // ones (resize), rather than leaving it to the compiler's heuristics, whose
    // choices flip with small changes in these large functions.
    template <bool orient_ = false, typename F, typename V>
    std::size_t insert_or_update(const char* const* seqs, const std::size_t* lens, std::size_t c, const F& f, const V& val, const Token& token);

public:

    // Constructs a streaming k-mer hash table to support up-to `max_sz` k-mers,
    // with a maximum load-factor of `lf`. Its bucket count is the fewest for
    // `max_sz` k-mers at that load, rounded up to a power of 2 iff `pow2_cap`.
    // Resizes use `resize_worker_c` threads.
    Streaming_Kmer_Hash_Table(std::size_t max_sz = 1 << 22, uint64_t resize_worker_c = 16, double lf = lf_default_, bool pow2_cap = false);

    // Returns whether a table constructed for `max_sz` k-mers at load factor
    // `lf`, rounded up to a power of 2 iff `pow2_cap`, is within the most
    // buckets supported.
    static bool fits(std::size_t max_sz, double lf, bool pow2_cap);

    Streaming_Kmer_Hash_Table(const Streaming_Kmer_Hash_Table&) = delete;
    Streaming_Kmer_Hash_Table(Streaming_Kmer_Hash_Table&&) = delete;
    Streaming_Kmer_Hash_Table& operator=(const Streaming_Kmer_Hash_Table&) = delete;
    Streaming_Kmer_Hash_Table& operator=(Streaming_Kmer_Hash_Table&&) = delete;

    ~Streaming_Kmer_Hash_Table();

    // Returns the capacity (in terms of keys) of the table.
    auto capacity() const { return cap_ * B; }

    // Returns the number of bytes of the table's buckets.
    std::size_t bucket_bytes() const { return cap_ * Store::bucket_bytes; }

    // Returns the size of the table.
    std::size_t size() const;

    // Returns the size of the overflow table.
    std::size_t overflow_size() const { return ov->size();  }

    // Inserts the k-mers of the sequence `seq` of length `len` into the set,
    // skipping the ones already present. Non-ACGT characters split the
    // sequence. Returns the number of k-mers inserted. The invoking user
    // thread must provide its token `token`.
    std::size_t insert(const char* seq, std::size_t len, const Token& token) requires (is_set_) { return insert_or_update(&seq, &len, 1, [](){}, 0, token); }

    // Inserts the k-mers of the sequences `seqs[0..c)` of lengths `lens[0..c)`
    // into the set, as `insert` does for one sequence. Batching sequences,
    // such as reads, lets their k-mers be described at once.
    std::size_t insert(const char* const* seqs, const std::size_t* lens, std::size_t c, const Token& token) requires (is_set_) { return insert_or_update(seqs, lens, c, [](){}, 0, token); }

    // Updates the associated value `v` of each k-mer of the sequence `seq` of
    // length `len` present in the map to `f(v)`, and inserts each absent one
    // with associated value `val`, one k-mer at a time, in an unspecified
    // order. Non-ACGT characters split the sequence. Returns the number of
    // k-mers inserted. The invoking user thread must provide its token
    // `token`.
    template <typename F, typename V = T_> std::size_t upsert(const char* seq, std::size_t len, const F& f, const V& val, const Token& token) requires (is_map_) { return insert_or_update(&seq, &len, 1, f, val, token); }

    // Upserts the k-mers of the sequences `seqs[0..c)` of lengths `lens[0..c)`
    // into the map, as `upsert` does for one sequence. Batching sequences,
    // such as reads, lets their k-mers be described at once.
    template <typename F, typename V = T_> std::size_t upsert(const char* const* seqs, const std::size_t* lens, std::size_t c, const F& f, const V& val, const Token& token) requires (is_map_) { return insert_or_update(seqs, lens, c, f, val, token); }

    // Upserts the k-mers of the sequences `seqs[0..c)` of lengths `lens[0..c)`
    // into the map with their orientations: updates the associated value `v`
    // of each k-mer to `f(v, fwd)`, from `T_()` for an absent one, one k-mer
    // at a time, in an unspecified order, where `fwd` is whether the k-mer
    // occurs in its canonical form, i.e. as is rather than as its reverse
    // complement. `f` must have no side effects. Non-ACGT characters split the
    // sequences. Returns the number of k-mers inserted. The invoking user
    // thread must provide its token `token`.
    template <typename F> std::size_t upsert_oriented(const char* const* seqs, const std::size_t* lens, std::size_t c, const F& f, const Token& token) requires (is_map_ && std::is_invocable_r_v<T_, const F&, T_, bool>) { return insert_or_update<true>(seqs, lens, c, f, T_(), token); }

    // Updates the associated value `v` of each k-mer of the sequence `seq` of
    // length `len` present in the map to `f(v)`, one k-mer at a time, in an
    // unspecified order, and skips the absent ones: it inserts none. Non-ACGT
    // characters split the sequence. In a concurrent map, each value updates
    // atomically, and `f` may be applied to it more than once, so it must have
    // no side effects. Returns the number of k-mers present. It takes no
    // locks, nor a token, and must not overlap operations on the map other
    // than `update`s.
    template <typename F> std::size_t update(const char* seq, std::size_t len, const F& f) requires (is_map_) { return update(&seq, &len, 1, f); }

    // Updates the values of the k-mers of the sequences `seqs[0..c)` of lengths
    // `lens[0..c)` present in the map, as `update` does for one sequence.
    // Batching sequences, such as reads, lets their k-mers be described at
    // once.
    template <typename F> std::size_t update(const char* const* seqs, const std::size_t* lens, std::size_t c, const F& f) requires (is_map_);

    // Looks up the k-mers of the sequence `seq` of length `len`, calling
    // `g(i, r)` for the k-mer at index `i`, with `r` as `true` in a set, or
    // its value in a map, iff it is present, in an unspecified order.
    // Non-ACGT characters split the sequence, and their k-mers are skipped.
    // Must not overlap updates that can resize the table.
    template <typename G> void find(const char* seq, std::size_t len, const G& g) const { find(&seq, &len, 1, [&](std::size_t, const std::size_t i, const find_ret_t& r){ g(i, r); }); }

    // Looks up the k-mers of the sequences `seqs[0..c)` of lengths
    // `lens[0..c)`, calling `g(q, i, r)` for the k-mer at index `i` of the
    // `q`'th sequence, as `find` does for one sequence. Iff `g` takes it,
    // calls `g(q, i, r, fwd)` instead, with `fwd` as whether the k-mer occurs
    // in its canonical form, as in `upsert_oriented`.
    template <typename G> void find(const char* const* seqs, const std::size_t* lens, std::size_t c, const G& g) const;

    // Looks up the k-mer of the `k` characters at `kmer`, in either
    // orientation. Returns its value in a map, or `true` in a set, iff it is
    // present; `std::nullopt`, or `false`, otherwise, as for a k-mer with a
    // non-ACGT character. Streams of k-mers go faster through the sequence
    // `find`s. Must not overlap updates that can resize the table.
    find_ret_t find(const char* kmer) const;

    // Returns the number of k-mers of the sequence `seq` of length `len`
    // present in the table. Must not overlap updates that can resize the
    // table.
    std::size_t count(const char* seq, std::size_t len) const { return count(&seq, &len, 1); }

    // Returns the number of k-mers of the sequences `seqs[0..c)` of lengths
    // `lens[0..c)` present in the table. Must not overlap updates that can
    // resize the table.
    std::size_t count(const char* const* seqs, const std::size_t* lens, std::size_t c) const;

    // Calls `f(x)` for each k-mer `x` in the set, or `f(x, v)` for each k-mer
    // `x` with its value `v` in the map, that lies in the share of worker
    // `w_id` of a group of `w_c` workers, in an unspecified order. The group's
    // shares partition the table, so that its workers, concurrent or not,
    // visit each k-mer once. Must not overlap operations that modify the
    // table.
    template <typename F> void for_each(const F& f, std::size_t w_c, std::size_t w_id) const;

    // Registers a new user and returns a unique token for it.
    Token register_user() { return Token(registered_thread_count++); }

    // Returns the percentages of keys in the table that are placed into their
    // primary bucket, their secondary buckets, and the overflow table. Their
    // minimizers are recomputed, by as many threads as resizes take.
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
inline Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Streaming_Kmer_Hash_Table(const std::size_t max_sz, const uint64_t resize_worker_c, const double lf, const bool pow2_cap):
      lf(lf)
    , cap_(cap_for(max_sz, lf, pow2_cap))
    , resize_th(capacity() * lf)
    , ov(new ov_t(std::max<std::size_t>(ov_min_cap, capacity() * of_default), resize_worker_c))
    , ov_new(nullptr)
    , registered_thread_count(0)
    , approx_sz(0)
    , checkpoint(max_user, Padded<uint64_t>(0))
    , checkpoint_step(step_for(resize_th))
    , resize_worker_c(resize_worker_c)
    , ov_sz(0)
    , ov_checkpoint(max_user, Padded<uint64_t>(0))
{
    if(cap_ > max_cap)
    {
        std::cerr << "A table of " << cap_ << " buckets exceeds the most supported, " << max_cap << ". Aborting.\n";
        std::abort();
    }

    set_ov_th();

    // Map the table up front, in parallel, rather than on first touches by the updates. The metadata comes zeroed:
    // checksum `0` denotes absence of a key.
    T.allocate(cap_, resize_worker_c);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::fits(const std::size_t max_sz, const double lf, const bool pow2_cap)
{
    // NB: The first test keeps the bucket count in range of `cap_for`'s conversion.
    return max_sz / lf / B <= max_cap && cap_for(max_sz, lf, pow2_cap) <= max_cap;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline Streaming_Kmer_Hash_Table<k, mt_, T_, l>::~Streaming_Kmer_Hash_Table()
{
    T.deallocate(cap_), delete(ov);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Store::allocate(const std::size_t n, const std::size_t t_c)
{
    T = allocate_huge<Bucket>(n);
    parallel_prefault(T, n * sizeof(Bucket), t_c);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Store::grow(const std::size_t n, const std::size_t m)
{
    T = grow_huge(T, n, m);
    return T != nullptr;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Store::deallocate(const std::size_t n)
{
    deallocate_huge(T, n);
    T = nullptr;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline unsigned char* Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Store::at(const std::size_t i, const std::size_t j) const
{
    return T[i].slot + j * slot_sz;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline Kmer<k> Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Store::key(const std::size_t i, const std::size_t j) const
{
    if constexpr(!packed)
        return Streaming_Kmer_Hash_Table::key(T[i].slot[j]);
    else
    {
        Kmer<k> x;
        std::memcpy(static_cast<void*>(&x), at(i, j), key_sz);
        return x;
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Store::val(const std::size_t i, const std::size_t j) const -> mapped_t
{
    if constexpr(is_set_)
        return (void)i, (void)j, mapped_t();
    else if constexpr(!packed)
        return T[i].slot[j].second;
    else
    {
        mapped_t v;
        std::memcpy(static_cast<void*>(&v), at(i, j) + key_sz, val_sz);
        return v;
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Store::set_val(const std::size_t i, const std::size_t j, const mapped_t& v) const
{
    if constexpr(is_set_)
        (void)i, (void)j, (void)v;
    else if constexpr(!packed)
        T[i].slot[j].second = v;
    else
        std::memcpy(at(i, j) + key_sz, static_cast<const void*>(&v), val_sz);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename F>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Store::update_val(const std::size_t i, const std::size_t j, const F& f) const
{
    if constexpr(!packed)
        update_in_place(T[i].slot[j].second, f);
    else if constexpr(!mt_)
        set_val(i, j, f(val(i, j)));
    else
    {
        // A packed value sits at offset `sizeof(Metadata) + j * slot_sz + key_sz` of its aligned bucket.
        static_assert((sizeof(Metadata) + key_sz) % std::atomic_ref<mapped_t>::required_alignment == 0 &&
                      slot_sz % std::atomic_ref<mapped_t>::required_alignment == 0, "Packed values are aligned for atomic updates.");
        update_in_place(*reinterpret_cast<mapped_t*>(at(i, j) + key_sz), f);
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Store::get(const std::size_t i, const std::size_t j) const -> flat_t
{
    if constexpr(!packed)
        return T[i].slot[j];
    else
        return flat_t(key(i, j), val(i, j));
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Store::put(const std::size_t i, const std::size_t j, const flat_t& e) const
{
    if constexpr(!packed)
        T[i].slot[j] = e;
    else
    {
        std::memcpy(at(i, j), static_cast<const void*>(&e.first), key_sz);
        std::memcpy(at(i, j) + key_sz, static_cast<const void*>(&e.second), val_sz);
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool write_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Store::prefetch_slot(const std::size_t i, const std::size_t j) const
{
    if constexpr(!packed)
        __builtin_prefetch(&T[i].slot[j], write_, 3);
    else
    {
        __builtin_prefetch(at(i, j), write_, 3);
        if constexpr(line % slot_sz != 0)  // The slot may cross into the next line.
            __builtin_prefetch(at(i, j) + slot_sz - 1, write_, 3);
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Store::prefetch_run(const std::size_t i, const std::size_t e) const
{
    if constexpr(!packed)
    {
        __builtin_prefetch(&T[i].slot[e], 1, 3);
        if(e + line_slots < B)
            __builtin_prefetch(&T[i].slot[e + line_slots], 1, 3);
    }
    else
    {
        __builtin_prefetch(at(i, e), 1, 3);
        if(e < run_lim)
            __builtin_prefetch(at(i, e) + line, 1, 3);
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Store::prefetch_frontier(const std::size_t i, const std::size_t e) const
{
    if constexpr(!packed)
    {
        if(e % line_slots == 0 && e + line_slots < B)
            __builtin_prefetch(&T[i].slot[e + line_slots], 1, 3);
    }
    else if((line_starts >> e) & 1)
        __builtin_prefetch(at(i, e) + line, 1, 3);
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
inline std::size_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::cap_for(const std::size_t max_sz, const double lf, const bool pow2)
{
    auto cap = std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(max_sz / lf / B)));
    while(static_cast<uint64_t>(cap * B * lf) <= max_sz)    // As `resize_th` takes it; it can round onto `max_sz`.
        cap++;

    return pow2 ? ceil_pow_2(cap) : cap;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline std::size_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::bucket_of(const uint64_t h, const std::size_t cap)
{
    assert(cap <= max_cap);
    return ((h & 0xFFFF'FFFFull) * cap) >> 32;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline std::size_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::bucket_of(const uint64_t h) const
{
    return bucket_of(h, cap_);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline std::size_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::bucket_of(const Kmer_Desc& d)
{
    return d.b;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::buckets_of(const uint64_t* const h, const std::size_t n, const std::size_t cap, std::size_t* const b)
{
    // NB: As widening 32-bit products, which vectorize.
    for(std::size_t j = 0; j < n; ++j)
        b[j] = (static_cast<uint64_t>(static_cast<uint32_t>(h[j])) * static_cast<uint32_t>(cap)) >> 32;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::buckets_of(const typename Kmer_Lanes<k, l>::vec h, const std::size_t cap) -> typename Kmer_Lanes<k, l>::vec
{
    typedef typename Kmer_Lanes<k, l>::vec vec;
    // NB: A single 32-bit multiply per k-mer: denser vector multiplies drop Cascade Lake cores to a lower AVX2 clock.
    const auto p = _mm256_mul_epu32(reinterpret_cast<__m256i>(h), _mm256_set1_epi64x(static_cast<long long>(cap)));
    return reinterpret_cast<vec>(_mm256_srli_epi64(p, 32));
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline std::size_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::pick_bucket(const uint64_t x, const std::size_t cap)
{
    return static_cast<std::size_t>((static_cast<unsigned __int128>(x) * cap) >> 64);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline std::size_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::filter_pos(const uint64_t h)
{
    const auto x = (h >> 32) & 0xFF;
    return x + (x == 0);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::may_spill(const Metadata& md, const uint64_t h)
{
    const auto i = filter_pos(h);
    return md.flags[i / 8] & (1u << (i % 8));
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::mark_spill(Metadata& md, const uint64_t h)
{
    // NB: Lock waiters CAS all of `flags[0]`, so marks go in under the bucket's lock where the bucket is shared.
    const auto i = filter_pos(h);
    md.flags[i / 8] |= (1u << (i % 8));
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::clear_spills(Metadata& md)
{
    md.flags[0] &= Metadata::lock_mask;
    std::fill(md.flags + 1, std::end(md.flags), 0);
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
    const auto n = bucket_size(b);
    uint64_t h[B];
    bucket_hashes(b, n, h);

    std::size_t p = 0, s = 0;
    for(std::size_t j = 0; j < n; ++j)
        p += (bucket_of(h[j]) == b), s += (bucket_of(h[j]) != b);

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
                // A waiter spins on loads, pausing, and tries the CAS only once the bucket looks unlocked: spinning
                // CASes would take the metadata line from the holder, which works on it; builds collapsed on hot
                // buckets at 56–128 threads.
                auto& w = (resize_ ? md_new(b) : md(b)).flags[0];
                const auto cur = __atomic_load_n(&w, __ATOMIC_RELAXED);
                if(!(cur & Metadata::lock_mask) && __sync_bool_compare_and_swap(&w, cur, cur | Metadata::lock_mask))
                    break;

                _mm_pause();
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
                assert(md(b).flags[0] & Metadata::lock_mask);
                md(b).flags[0] &= ~Metadata::lock_mask;
            }
            else
            {
                assert(md_new(b).flags[0] & Metadata::lock_mask);
                md_new(b).flags[0] &= ~Metadata::lock_mask;
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
        return !p ? null_val : val_t(*p);
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename F>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::update_in_place(mapped_t& v, const F& f)
{
    if constexpr(mt_)
    {
        static_assert(std::atomic_ref<mapped_t>::is_always_lock_free, "Values of a concurrent map update in place without locks.");
        assert(reinterpret_cast<uintptr_t>(&v) % std::atomic_ref<mapped_t>::required_alignment == 0);

        const std::atomic_ref<mapped_t> a(v);
        auto x = a.load(std::memory_order_relaxed);
        while(!a.compare_exchange_weak(x, f(x), std::memory_order_relaxed))
            ;
    }
    else
        v = f(v);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename F>
bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::update_in_overflow(const Kmer<k> key, const F& f)
{
    const auto p = ov->find(key);
    if(p)
        update_in_place(*p, f);

    return p != nullptr;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::place_at(const flat_t& e, uint8_t c, std::size_t i, std::size_t j)
{
    T.put(i, j, e), md(i).cs[j] = c;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::place_at_resize(const flat_t& e, const uint8_t c, const std::size_t i, const std::size_t j)
{
    T_new.put(i, j, e), md_new(i).cs[j] = c;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::try_insert_at_resize(const flat_t& e, const uint8_t c, const uint64_t h, const std::size_t i)
{
    lock<true>(i);

    const auto c_vec = load_cs(md_new(i));
    const auto zero = zero_cs();
    const auto cmp = eq_mask(c_vec, zero);

    const auto j = std::countr_zero(cmp);
    if(j == B)
    {
        mark_spill(md_new(i), h);
        unlock<true>(i);
        return false;
    }

    place_at_resize(e, c, i, j);
    unlock<true>(i);
    return true;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::insert_at_resize(const flat_t& x, const uint8_t c, const uint64_t h, const bool sec, std::size_t p, std::size_t q)
{
    const auto zero = zero_cs();
    const auto full = [&](const std::size_t i){ return eq_mask(load_cs(md_new(i)), zero) == 0; };

    // A full primary bucket's spill mark is tested without its lock, and set with it.
    const auto b = bucket_of(h);
    if(!full(b))
    {
        if(try_insert_at_resize(x, c, h, b))
            return false;
    }
    else if(!may_spill(md_new(b), h))
        lock<true>(b), mark_spill(md_new(b), h), unlock<true>(b);

    if(!sec)
        std::tie(p, q) = secondaries(h);
    if(p > q)
        std::swap(p, q);

    const bool locked = !(full(p) && full(q));
    if(locked)
        lock_ordered<true>(p, q);

    const auto c_vec_p = load_cs(md_new(p));
    const auto c_vec_q = load_cs(md_new(q));

    const auto cmp_p = eq_mask(c_vec_p, zero);
    const auto cmp_q = eq_mask(c_vec_q, zero);

    const auto sz_p = static_cast<std::size_t>(std::countr_zero(cmp_p));
    const auto sz_q = static_cast<std::size_t>(std::countr_zero(cmp_q));
    assert(std::countr_zero(cmp_p) + std::countl_one(cmp_p) == B), assert(std::countr_zero(cmp_q) + std::countl_one(cmp_q) == B);

    if(sz_p < sz_q)
    {
        place_at_resize(x, c, p, sz_p);
        unlock_ordered<true>(p, q);
    }
    else if(sz_q < B)
    {
        place_at_resize(x, c, q, sz_q);
        unlock_ordered<true>(p, q);
    }
    else
    {
        if(locked)
            unlock_ordered<true>(p, q);

        bool ov_full;
        if constexpr(is_set_)
            ov_new->insert(x, Overflow_Set_Val(), ov_full);
        else
            ov_new->insert(x.first, x.second, ov_full);

        if(ov_full)
        {
            std::cerr << "Overflow table of capacity " << ov_new->capacity() << " is full in a resize. Aborting.\n";
            std::abort();
        }

        return true;
    }

    return false;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::add_ov_key(const uint16_t id)
{
    auto& cp = ov_checkpoint[id].unwrap();
    if(KACHE_HASH_LIKELY(++cp < ov_step.load(std::memory_order_relaxed)))
        return false;

    const auto sz = (ov_sz += cp);
    cp = 0;
    return sz >= ov_th.load(std::memory_order_relaxed);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline std::size_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::new_ov_size() const
{
    const auto base = std::max<std::size_t>(ov_min_cap, capacity() * of_default);
    return std::max<std::size_t>(base, ov->capacity() * ov_t::lf_default());
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::reset_ov_size(const uint64_t sz)
{
    ov_sz = sz;
    std::for_each(ov_checkpoint.begin(), ov_checkpoint.end(), [](auto& c){ c = 0; });
    set_ov_th();
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::set_ov_th()
{
    ov_th = static_cast<uint64_t>(ov->capacity() * ov_lf);
    ov_step = std::clamp<uint64_t>(ov_th / (8 * max_user), 1, 1024);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::try_grow_ov(const ov_t* const ov_seen, const bool full)
{
    if constexpr(mt_)   table_lock.lock();

    if(ov == ov_seen && (full || ov_sz >= ov_th))
    {
        // Re-insert the keys into a table of twice the capacity, in parallel.
        auto* ov_2 = new ov_t(ov->capacity() * 2 * ov_t::lf_default(), resize_worker_c);
        std::vector<std::thread> W;
        W.reserve(resize_worker_c);
        for(uint64_t w = 0; w < resize_worker_c; ++w)
            W.emplace_back([&](const uint64_t w_id)
            {
                auto it = ov->iterator(resize_worker_c, w_id);
                Kmer<k> key;
                ov_val_t val;
                bool full_2;
                while(it.next(key, val))
                    ov_2->insert(key, val, full_2);

                ov->release(resize_worker_c, w_id);
            }, w);

        std::for_each(W.begin(), W.end(), [](auto& w){ w.join(); });
        std::swap(ov, ov_2);
        delete(ov_2);
        set_ov_th();
    }

    if constexpr(mt_)   table_lock.unlock();
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
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::bucket_hashes(const std::size_t b, const std::size_t n, uint64_t* const h) const
{
    // NB: Returning on empty buckets also keeps GCC from warning that `x` may be read uninitialized.
    if(n == 0)
        return;

    uint64_t x[B];
    for(std::size_t j = 0; j < n; ++j)
        x[j] = T.key(b, j).as_int();

    Kmer_Lanes<k, l>::minimizer_hashes(x, n, h);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::resize()
{
    const auto old_cap = cap_;
    if(2 * old_cap > max_cap)
    {
        std::cerr << "Growing the table to " << 2 * old_cap << " buckets exceeds the most supported, " << max_cap << ". Aborting.\n";
        std::abort();
    }

    if(!T.grow(old_cap, 2 * old_cap))   // The new buckets come zeroed: checksum `0` denotes absence of a key.
    {
        std::cerr << "Growing the table to " << 2 * old_cap * B << " slots failed. Aborting.\n";
        std::abort();
    }

    cap_ *= 2;
    T_new = T;  // The resize-insertion helpers act on this.
    ov_new = new ov_t(new_ov_size(), resize_worker_c);
    std::atomic<uint64_t> ov_new_c{0};  // Number of keys moved into the new overflow table.

    // In the first phase, the keys at their primary buckets `b` move to `2b` or `2b + 1`, without locks; the keys at
    // secondary buckets are set aside, and with the overflow table's, re-inserted with locks in the second phase.
    // `2b` and `2b + 1` are to be emptied first, so the old buckets go in rounds from the top, over halving ranges:
    // a round's `[s, e)`, with `s = ceil(e / 2)`, moves to `[2s, 2e)`, which the earlier rounds emptied, or are new.
    // The workers split each round's range. Bucket 0 goes last, keeping its keys for `2b = 0` compacted.
    const auto worker_c = resize_worker_c;
    struct Held { flat_t x; uint64_t h; uint8_t c; };   // A key set aside, with its minimizer hash.
    std::vector<std::vector<Held>> held_of(worker_c);   // The keys each worker set aside.
    std::vector<std::size_t> held_pre(worker_c + 1);    // Prefix sums of their numbers, once all are set aside.
    std::atomic<std::size_t> held_next{0};  // Next chunk of them to re-insert.
    const auto seal = [&]() noexcept
    {
        for(std::size_t w = 0; w < worker_c; ++w)
            held_pre[w + 1] = held_pre[w] + held_of[w].size();
    };
    std::barrier phase_end(static_cast<std::ptrdiff_t>(worker_c), seal);
    std::barrier round_end(static_cast<std::ptrdiff_t>(worker_c));  // Ends the first phase's rounds.
    std::vector<std::thread> W;
    W.reserve(worker_c);
    for(uint64_t i = 0; i < worker_c; ++i)
        W.emplace_back([&](const auto w_id)
        {
            std::vector<Held> held; // This worker's keys at secondary buckets: local while filled, as workers pushing
                                    // into vectors whose headers share cache lines would contend for them.

            // Room for a quarter of the slots of the worker's share of the old buckets, typically more than its keys
            // set aside: only what's touched gets backed, and the vector need not copy itself to grow.
            held.reserve((old_cap / worker_c + 1) * B / 4);

            // Moves the keys of the old bucket `b` at their primary buckets to `2b` and `2b + 1`, compacting the ones
            // at `2b` if it is `b` itself; sets the rest aside.
            uint64_t h_b[B];    // Minimizer hashes of a bucket's keys.
            const auto split = [&](const std::size_t b)
            {
                const auto lb = 2 * b, hb = lb + 1;
                std::size_t lo = 0, hi = 0; // Numbers of keys moved to `lb`, and to `hb`.
                const auto n = bucket_size(b);
                bucket_hashes(b, n, h_b);
                for(std::size_t j = 0; j < n; ++j)
                {
                    const uint8_t c = md(b).cs[j];
                    const auto h = h_b[j];
                    if(bucket_of(h, old_cap) != b)
                    {
                        held.push_back({T.get(b, j), h, c});
                        continue;
                    }

                    if(bucket_of(h) == hb)
                        T.put(hb, hi, T.get(b, j)), md(hb).cs[hi] = c, hi++;
                    else
                    {
                        if(lb != b || lo != j)
                            T.put(lb, lo, T.get(b, j)), md(lb).cs[lo] = c;
                        lo++;
                    }
                }

                for(std::size_t j = (lb == b ? lo : 0); j < n; ++j)
                    md(b).cs[j] = 0;

                clear_spills(md(b));    // The second phase marks the spills anew.
            };

            for(std::size_t e = old_cap; e > 1; )
            {
                const auto s = (e + 1) / 2;
                const auto b_s = s + (e - s) * w_id / worker_c, b_e = s + (e - s) * (w_id + 1) / worker_c;
                for(std::size_t b = b_s; b < b_e; ++b)
                    split(b);

                round_end.arrive_and_wait();
                e = s;
            }

            if(w_id == 0)
                split(0);

            held_of[w_id] = std::move(held);
            phase_end.arrive_and_wait();

            // The rest go through a pipeline of `pipe_sz` keys: a key's new primary bucket is prefetched as it enters;
            // `stage_d` keys later, its checksums are read to prefetch the slot it likely takes, or if it is full, its
            // secondary buckets and its overflow slot; the key is inserted as it leaves.
            constexpr std::size_t pipe_sz = 32, stage_d = 16;
            struct Pending { flat_t x; uint64_t h; std::size_t p, q; uint8_t c; bool sec; };
            Pending pipe[pipe_sz];
            std::size_t pipe_c = 0;
            uint64_t ov_n = 0;  // Number of keys this worker moves into the new overflow table.
            const auto stage = [&](Pending& e)
            {
                const auto i = bucket_of(e.h);
                const auto empty = eq_mask(load_cs(md(i)), zero_cs());
                if(empty)
                    T.template prefetch_slot<true>(i, std::countr_zero(empty));
                else
                {
                    std::tie(e.p, e.q) = secondaries(e.h), e.sec = true;
                    __builtin_prefetch(&md(e.p), 1, 3), __builtin_prefetch(&md(e.q), 1, 3);
                    ov_new->template prefetch<true>(key(e.x));
                }
            };
            const auto insert = [&](const Pending& e){ ov_n += insert_at_resize(e.x, e.c, e.h, e.sec, e.p, e.q); };
            const auto push = [&](const flat_t& x, const uint64_t h, const uint8_t c)
            {
                const auto& meta = md(bucket_of(h));
                __builtin_prefetch(meta.cs, 1, 3);
                if constexpr(sizeof(Metadata) > L1_CACHE_LINE_SIZE)
                    __builtin_prefetch(meta.flags, 1, 3);

                if(pipe_c >= stage_d)
                    stage(pipe[(pipe_c - stage_d) % pipe_sz]);
                auto& e = pipe[pipe_c % pipe_sz];
                if(pipe_c >= pipe_sz)
                    insert(e);
                e = Pending{x, h, 0, 0, c, false};
                pipe_c++;
            };

            // The keys set aside go in chunks, taken dynamically from all the workers'.
            constexpr std::size_t chunk_sz = 16384;
            const auto held_c = held_pre[worker_c];
            for(std::size_t ch; (ch = held_next++) * chunk_sz < held_c; )
            {
                // The `x`'th key set aside is the `v`'th worker's.
                const auto s = ch * chunk_sz, e = std::min(held_c, s + chunk_sz);
                std::size_t v = std::upper_bound(held_pre.cbegin(), held_pre.cend(), s) - held_pre.cbegin() - 1;
                for(std::size_t x = s; x < e; ++x)
                {
                    while(x >= held_pre[v + 1])
                        v++;

                    const auto& y = held_of[v][x - held_pre[v]];
                    push(y.x, y.h, y.c);
                }
            }

            // The overflow table's keys go in batches, whose minimizer hashes are recomputed at once, and their
            // checksums one by one.
            auto it = ov->iterator(worker_c, w_id);
            constexpr std::size_t batch_sz = 32;
            Kmer<k> o_key[batch_sz];
            ov_val_t o_val[batch_sz];
            uint64_t o_x[batch_sz], o_h[batch_sz];
            std::size_t o_c = 0;    // Number of keys in the batch.
            const auto push_batch = [&]()
            {
                for(std::size_t j = 0; j < o_c; ++j)
                    o_x[j] = o_key[j].as_int();
                Kmer_Lanes<k, l>::minimizer_hashes(o_x, o_c, o_h);

                for(std::size_t j = 0; j < o_c; ++j)
                {
                    if constexpr(is_set_)
                        push(o_key[j], o_h[j], kmer_checksum(o_x[j]));
                    else
                        push(flat_t(o_key[j], o_val[j]), o_h[j], kmer_checksum(o_x[j]));
                }

                o_c = 0;
            };

            while(it.next(o_key[o_c], o_val[o_c]))
                if(++o_c == batch_sz)
                    push_batch();
            push_batch();

            // The pipeline's last keys.
            for(std::size_t e = (pipe_c > stage_d ? pipe_c - stage_d : 0); e < pipe_c; ++e)
                stage(pipe[e % pipe_sz]);
            for(std::size_t e = (pipe_c > pipe_sz ? pipe_c - pipe_sz : 0); e < pipe_c; ++e)
                insert(pipe[e % pipe_sz]);

            ov_new_c += ov_n;

            ov->release(worker_c, w_id);
        }, i);

    std::for_each(W.begin(), W.end(), [](auto& w){ w.join(); });

    std::swap(ov, ov_new);
    delete(ov_new);
    T_new = Store(), ov_new = nullptr;
    reset_ov_size(ov_new_c);

    resize_th = capacity() * lf;
    checkpoint_step = step_for(resize_th);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool lanes_, bool cont_, bool buckets_, typename Pre_Meta_, typename Pre_Slot_, typename Op_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::stream(const char* const* const seqs, const std::size_t* const lens, const std::size_t c, const std::size_t cap, const Pre_Meta_& pre_meta, const Pre_Slot_& pre_slot, const Op_& op)
{
    if constexpr(!lanes_)    // All stretches through one pipeline, in groups.
    {
        constexpr std::size_t group_sz = 64;
        Stretch S[group_sz];
        std::size_t s_c = 0;
        for(std::size_t q = 0; q < c; ++q)
        {
            const char* const seq = seqs[q];
            const auto len = lens[q];
            for(std::size_t s = 0; s + k <= len; )  // `s` is the start of the current stretch of ACGT bases.
            {
                const auto e = s + DNA_Utility::acgt_prefix(seq + s, len - s);   // End of the stretch.

                if(e - s >= k)
                {
                    S[s_c++] = {seq + s, e - s, q, s};
                    if(s_c == group_sz)
                        stream_stretches<buckets_>(S, s_c, cap, pre_meta, pre_slot, op), s_c = 0;
                }

                s = e + 1;
            }
        }

        if(s_c > 0)
            stream_stretches<buckets_>(S, s_c, cap, pre_meta, pre_slot, op);
    }
    else
        stream_lanes_all<cont_, buckets_>(seqs, lens, c, cap, pre_meta, pre_slot, op);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool cont_, bool buckets_, typename Pre_Meta_, typename Pre_Slot_, typename Op_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::stream_lanes_all(const char* const* const seqs, const std::size_t* const lens, const std::size_t c, const std::size_t cap, const Pre_Meta_& pre_meta, const Pre_Slot_& pre_slot, const Op_& op)
{
    // Groups of stretches gathered to go through the lanes in one pipeline: a long stretch split over the lanes on
    // its own, or short ones a lane each.
    Stretch G[lane_group_c * lane_c];
    std::size_t G_c[lane_group_c];  // Number of stretches in each group.
    std::size_t g_c = 0;    // Number of groups gathered, the last one possibly partial.
    const auto flush = [&]()
    {
        stream_lanes<cont_, buckets_>(G, G_c, g_c, cap, pre_meta, pre_slot, op);
        g_c = 0;
    };
    const auto open_group = [&]()   // Makes the last group a fresh one, flushing the groups first if full.
    {
        if(g_c == lane_group_c)
            flush();
        G_c[g_c++] = 0;
    };

    for(std::size_t q = 0; q < c; ++q)
    {
        const char* const seq = seqs[q];
        const auto len = lens[q];
        for(std::size_t s = 0; s + k <= len; )  // `s` is the start of the current stretch of ACGT bases.
        {
            const auto e = s + DNA_Utility::acgt_prefix(seq + s, len - s);   // End of the stretch.

            if(e - s >= k)
            {
                const auto n = e - s - (k - 1);
                if(n >= lane_min_kmers)    // Split over the lanes on its own, in a group placed before any partial one.
                {
                    if(g_c > 0 && G_c[g_c - 1] < lane_c)
                    {
                        Stretch P[lane_c];
                        const auto p_c = G_c[g_c - 1];
                        std::copy(G + (g_c - 1) * lane_c, G + (g_c - 1) * lane_c + p_c, P);
                        g_c--;
                        open_group();
                        for(std::size_t j = 0; j < lane_c; ++j)
                        {
                            const auto b = n * j / lane_c, b_e = n * (j + 1) / lane_c;
                            G[(g_c - 1) * lane_c + j] = {seq + s + b, b_e - b + (k - 1), q, s + b};
                        }
                        G_c[g_c - 1] = lane_c;
                        open_group();
                        std::copy(P, P + p_c, G + (g_c - 1) * lane_c);
                        G_c[g_c - 1] = p_c;
                    }
                    else
                    {
                        open_group();
                        for(std::size_t j = 0; j < lane_c; ++j)
                        {
                            const auto b = n * j / lane_c, b_e = n * (j + 1) / lane_c;
                            G[(g_c - 1) * lane_c + j] = {seq + s + b, b_e - b + (k - 1), q, s + b};
                        }
                        G_c[g_c - 1] = lane_c;
                    }
                }
                else
                {
                    if(g_c == 0 || G_c[g_c - 1] == lane_c)
                        open_group();
                    G[(g_c - 1) * lane_c + G_c[g_c - 1]++] = {seq + s, e - s, q, s};
                }
            }

            s = e + 1;
        }
    }

    // A last group of few stretches goes through the batch describer instead.
    std::size_t p_c = 0;
    if(g_c > 0 && G_c[g_c - 1] < lane_c / 2)
        p_c = G_c[--g_c];

    const auto* const P = G + g_c * lane_c;
    if(g_c > 0)
        flush();

    for(std::size_t j = 0; j < p_c; ++j)
        stream_batched<buckets_>(P[j].t, P[j].len - (k - 1), P[j].q, P[j].base, cap, pre_meta, pre_slot, op);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool buckets_, typename Pre_Meta_, typename Pre_Slot_, typename Op_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::stream_batched(const char* const t, const std::size_t n, const std::size_t q, const std::size_t base, const std::size_t cap, const Pre_Meta_& pre_meta, const Pre_Slot_& pre_slot, const Op_& op)
{
    typedef Kmer_Batch<k, l> batch_t;
    constexpr bool orient_ = oriented_op<Op_>;
    typedef std::conditional_t<orient_, typename batch_t::Oriented_Chunk, typename batch_t::Chunk> chunk_t;
    batch_t kb;
    chunk_t ch[2];  // The chunk being operated on, and the next; each is refilled two chunks ahead once done.
    std::size_t bk[2][buckets_ ? batch_t::C : 1];   // Iff `buckets_`, the primary buckets of the chunks' k-mers.

    const auto desc = [&](const std::size_t x, const std::size_t j){ const auto& c = ch[x]; return Kmer_Desc{Kmer<k>(c.key[j]), c.h[j], c.c[j], buckets_ ? bk[x][j] : 0}; };
    const auto fill = [&](const std::size_t x)
    {
        kb.next(ch[x]);
        if constexpr(buckets_)
            buckets_of(ch[x].h, ch[x].n, cap, bk[x]);
        for(std::size_t j = 0; j < ch[x].n; ++j)
            pre_meta(desc(x, j));
    };

    kb.reset(t, n);
    fill(0), fill(1);
    std::size_t first = base;   // Index of the current chunk's first k-mer.
    for(std::size_t cur = 0; ch[cur].n > 0; cur ^= 1)
    {
        const auto n_cur = ch[cur].n, n_next = ch[cur ^ 1].n;
        for(std::size_t j = 0; j < n_cur; ++j)
        {
            const auto x = j + slot_dist;
            if(x < n_cur)
                pre_slot([&](){ return desc(cur, x); });
            else if(x - n_cur < n_next)
                pre_slot([&](){ return desc(cur ^ 1, x - n_cur); });

            if constexpr(orient_)
                op(desc(cur, j), q, first + j, ch[cur].o[j]);
            else
                op(desc(cur, j), q, first + j);
        }

        first += n_cur;
        fill(cur);
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool buckets_, typename Pre_Meta_, typename Pre_Slot_, typename Op_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::stream_stretches(const Stretch* const S, const std::size_t c, const std::size_t cap, const Pre_Meta_& pre_meta, const Pre_Slot_& pre_slot, const Op_& op)
{
    typedef Kmer_Batch<k, l> batch_t;
    constexpr bool orient_ = oriented_op<Op_>;
    typedef std::conditional_t<orient_, typename batch_t::Oriented_Chunk, typename batch_t::Chunk> chunk_t;
    constexpr std::size_t depth = 3;    // Number of chunks in flight: the one operated on, and the next ones, whose metadata is being prefetched. Short stretches make short chunks, so two ahead keep the prefetches early enough.
    batch_t kb;
    chunk_t ch[depth];  // Each is refilled `depth` chunks ahead once done.
    std::size_t ch_s[depth];    // Stretch of each chunk; chunks don't span stretches.
    std::size_t ch_first[depth];    // Index of each chunk's first k-mer in its stretch.

    std::size_t s_next = 0; // Next stretch to describe.
    std::size_t s_cur = 0;  // Stretch being described.
    std::size_t at = 0; // Number of k-mers of the stretch described.
    bool active = false;    // Whether a stretch is being described.

    std::size_t bk[depth][buckets_ ? batch_t::C : 1];   // Iff `buckets_`, the primary buckets of the chunks' k-mers.

    const auto desc = [&](const std::size_t x, const std::size_t j){ const auto& c = ch[x]; return Kmer_Desc{Kmer<k>(c.key[j]), c.h[j], c.c[j], buckets_ ? bk[x][j] : 0}; };
    const auto fill = [&](const std::size_t x)
    {
        while(true)
        {
            if(!active)
            {
                if(s_next == c)
                {
                    ch[x].n = 0;
                    return;
                }

                s_cur = s_next++, at = 0, active = true;
                kb.reset(S[s_cur].t, S[s_cur].len - (k - 1));
            }

            if(kb.next(ch[x]) == 0)
            {
                active = false;
                continue;
            }

            ch_s[x] = s_cur, ch_first[x] = at, at += ch[x].n;
            if constexpr(buckets_)
                buckets_of(ch[x].h, ch[x].n, cap, bk[x]);
            for(std::size_t j = 0; j < ch[x].n; ++j)
                pre_meta(desc(x, j));
            return;
        }
    };

    for(std::size_t x = 0; x < depth; ++x)
        fill(x);
    for(std::size_t cur = 0; ch[cur].n > 0; cur = (cur + 1) % depth)
    {
        const auto nxt = (cur + 1) % depth;
        const auto n_cur = ch[cur].n, n_next = ch[nxt].n;
        const auto& st = S[ch_s[cur]];
        for(std::size_t j = 0; j < n_cur; ++j)
        {
            const auto x = j + slot_dist;
            if(x < n_cur)
                pre_slot([&](){ return desc(cur, x); });
            else if(x - n_cur < n_next)
                pre_slot([&](){ return desc(nxt, x - n_cur); });

            if constexpr(orient_)
                op(desc(cur, j), st.q, st.base + ch_first[cur] + j, ch[cur].o[j]);
            else
                op(desc(cur, j), st.q, st.base + ch_first[cur] + j);
        }

        fill(cur);
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool cont_, bool buckets_, typename Pre_Meta_, typename Pre_Slot_, typename Op_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::stream_lanes(const Stretch* const S, const std::size_t* const G_c, const std::size_t g_c, const std::size_t cap, const Pre_Meta_& pre_meta, const Pre_Slot_& pre_slot, const Op_& op)
{
    typedef Kmer_Lanes<k, l> lanes_t;
    constexpr auto T_c = lanes_t::T;
    constexpr bool orient_ = oriented_op<Op_>;
    static_assert(T_c * lane_c <= 256);
    static_assert(slot_dist >= T_c, "The slot prefetches a chunk ahead are to cover the chunk.");

    lanes_t kl;
    std::size_t g_d = 0;    // Group being described.
    const auto start = [&](const std::size_t g)
    {
        const char* t[lane_c];
        std::size_t len[lane_c];
        for(std::size_t j = 0; j < G_c[g]; ++j)
            t[j] = S[g * lane_c + j].t, len[j] = S[g * lane_c + j].len;

        kl.reset(t, len, G_c[g]);
    };

    start(0);

    // `lane_depth` chunks: one being operated on, and the next ones, whose metadata is prefetched once they are
    // described, before the operations on the one ahead of them. A chunk's stretches are its group's.
    typename lanes_t::Chunk ch[lane_depth];
    std::size_t ch_g[lane_depth];   // Group of each chunk.
    typename lanes_t::vec bk[lane_depth][buckets_ ? T_c : 1];  // Iff `buckets_`, the primary buckets of the chunks' k-mers.
    const auto desc = [&](const std::size_t x, const std::size_t j, const std::size_t t)
    {
        return Kmer_Desc{Kmer<k>(ch[x].key[t][j]), ch[x].h[t][j], static_cast<uint8_t>(ch[x].c[t][j]), buckets_ ? static_cast<std::size_t>(bk[x][t][j]) : 0};
    };
    [[maybe_unused]] const auto fwd_of = [&](const std::size_t x, const std::size_t j, const std::size_t t) -> bool { return (ch[x].c[t][j] >> 8) & 1; };
    const auto fill_and_prefetch = [&](const std::size_t x) -> bool
    {
        while(!describe_lanes<orient_>(kl, ch[x]))  // The group is done; describe the next one, if any.
        {
            if(g_d + 1 == g_c)
                return false;
            start(++g_d);
        }

        ch_g[x] = g_d;
        const auto c = G_c[g_d];
        if(ch[x].any)
        {
            if constexpr(buckets_)
                for(std::size_t t = 0; t < T_c; ++t)
                    bk[x][t] = buckets_of(ch[x].h[t], cap);

            for(std::size_t j = 0; j < c; ++j)
            {
                const auto [t_b, t_e] = ch[x].range(j);
                for(std::size_t t = t_b; t < t_e; ++t)
                    pre_meta(desc(x, j, t));
            }
        }

        return true;
    };

    // A chunk's slot prefetches, all at once a chunk ahead. With two chunks in flight, their metadata was prefetched
    // just before.
    const auto pre_slots = [&](const std::size_t x)
    {
        if constexpr(requires { pre_slot.active(); })
            if(!pre_slot.active())
                return;

        if(ch[x].any)
            for(std::size_t j = 0; j < G_c[ch_g[x]]; ++j)
            {
                const auto [t_b, t_e] = ch[x].range(j);
                for(std::size_t t = t_b; t < t_e; ++t)
                    pre_slot([&](){ return desc(x, j, t); });
            }
    };

    std::size_t done[lane_c] = {};  // With `cont_`, the number of each lane's first k-mers in the current chunk operated on with the chunk before, as the end of its last run.

    bool live[lane_depth];
    for(std::size_t x = 0; x < lane_depth; ++x)
        live[x] = fill_and_prefetch(x);
    pre_slots(0);
    for(std::size_t cur = 0; live[cur]; cur = (cur + 1 == lane_depth ? 0 : cur + 1))
    {
        const auto cur_next = (cur + 1 == lane_depth ? 0 : cur + 1);   // The next chunk.
        if(live[cur_next])
            pre_slots(cur_next);

        const auto g = ch_g[cur];
        const auto c = G_c[g];
        const Stretch* const S_g = S + g * lane_c;  // The chunk's stretches.
        if(ch[cur].any)
            for(std::size_t j = 0; j < c; ++j)
            {
                const auto [t_s, t_e] = ch[cur].range(j);
                const auto t_b = t_s + done[j];
                const auto base = S_g[j].base + static_cast<std::size_t>(ch[cur].first);
                done[j] = 0;

                for(std::size_t t = t_b; t < t_e; ++t)
                    if constexpr(orient_)
                        op(desc(cur, j, t), S_g[j].q, base + t, fwd_of(cur, j, t));
                    else
                        op(desc(cur, j, t), S_g[j].q, base + t);

                // Continue the lane's last run into its k-mers of the next chunk, if of the same group.
                if constexpr(cont_)
                    if(t_e > t_s && t_e == T_c && live[cur_next] && ch_g[cur_next] == g && ch[cur_next].any)
                    {
                        const auto [u_b, u_e] = ch[cur_next].range(j);
                        const auto base_next = S_g[j].base + static_cast<std::size_t>(ch[cur_next].first);
                        const auto h_last = ch[cur].h[t_e - 1][j];
                        std::size_t u = u_b;
                        for(; u < u_e && ch[cur_next].h[u][j] == h_last; ++u)
                            if constexpr(orient_)
                                op(desc(cur_next, j, u), S_g[j].q, base_next + u, fwd_of(cur_next, j, u));
                            else
                                op(desc(cur_next, j, u), S_g[j].q, base_next + u);

                        done[j] = u - u_b;
                    }
            }

        live[cur] = fill_and_prefetch(cur);
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool orient_>
bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::describe_lanes(Kmer_Lanes<k, l>& kl, typename Kmer_Lanes<k, l>::Chunk& out)
{
    return kl.template next<orient_>(out);
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
        if(T.key(i, j) == key)
            return j;
    }

    return B;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline uint32_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::probe(const Kmer<k> key, const uint8_t c, const std::size_t i, bool& full) const
{
    assert(i < cap_);

    auto cmp = eq_mask(load_cs(md(i)), splat_cs(c));
    while(cmp)
    {
        const auto j = std::countr_zero(cmp);
        cmp &= (cmp - 1);
        if(T.key(i, j) == key)
            return j;
    }

    full = (md(i).cs[B - 1] != 0);
    return B;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
std::pair<std::size_t, std::size_t> Streaming_Kmer_Hash_Table<k, mt_, T_, l>::secondaries(const uint64_t h) const
{
    return {pick_bucket(XXH3_64bits_withSeed(&h, sizeof(h), 0), cap_), pick_bucket(XXH3_64bits_withSeed(&h, sizeof(h), 1llu << 63), cap_)};
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool write_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::prefetch_meta(const Kmer_Desc& d) const
{
    const auto& meta = md(write_ ? bucket_of(d.h) : bucket_of(d));
    __builtin_prefetch(meta.cs, write_, 3);
    if constexpr(write_ && sizeof(Metadata) > L1_CACHE_LINE_SIZE)  // Updates lock through the flags too.
        __builtin_prefetch(meta.flags, write_, 3);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline std::pair<std::size_t, std::size_t> Streaming_Kmer_Hash_Table<k, mt_, T_, l>::secondaries(const uint64_t h, Sec_Cache& sc) const
{
    if(sc.cap != cap_ || sc.h != h)
    {
        sc.h = h, sc.cap = cap_;
        std::tie(sc.p, sc.q) = secondaries(h);
    }

    return {sc.p, sc.q};
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool write_>
void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::prefetch_spill(const Kmer<k> key, const uint64_t h, Sec_Cache& sc) const
{
    // The secondary buckets' metadata, once per minimizer run: with the overflow-table slot.
    if(sc.h != h || sc.cap != cap_)
    {
        const auto [p, q] = secondaries(h, sc);
        __builtin_prefetch(&md(p), write_, 3), __builtin_prefetch(&md(q), write_, 3);
    }

    ov->template prefetch<write_>(key);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool write_, bool insert_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::prefetch_slot(const Kmer_Desc& d, Sec_Cache& sc) const
{
    const auto i = (insert_ ? bucket_of(d.h) : bucket_of(d));
    const auto c_vec = load_cs(md(i));
    const auto cand = eq_mask(c_vec, splat_cs(d.c));
    if(!cand)
    {
        const auto empty = eq_mask(c_vec, zero_cs());
        if(!empty)  // The operation goes beyond the primary bucket.
        {
            if constexpr(!insert_)
                if(!may_spill(md(i), d.h))
                    return;

            return prefetch_spill<insert_>(d.key, d.h, sc);
        }

        // Updates predicted to insert prefetch the line after their predicted slot's too: a run's insertions cross
        // into it, and with the insertion-frontier prefetches alone, their stores still filled the store buffer.
        if constexpr(insert_)
            T.prefetch_run(i, std::countr_zero(empty));

        return;
    }

    T.template prefetch_slot<write_>(i, std::countr_zero(cand));
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool write_>
inline bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Lookup_Pre_Slot<write_>::active() const
{
    return hit;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool write_>
template <typename Make_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Lookup_Pre_Slot<write_>::operator()(const Make_& make) const
{
    if(hit)
        M.template prefetch_slot<write_, false>(make(), sc);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename F, typename V>
bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::ov_update(const Kmer<k> key, const F& f, const V& val, bool& full)
{
    if constexpr(is_set_)
        return (void)f, (void)val, !ov->insert(key, Overflow_Set_Val(), full);
    else
        return !ov->upsert(key, val, f, full);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool orient_, typename F, typename V>
__attribute__((flatten)) inline std::size_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::insert_or_update(const char* const* const seqs, const std::size_t* const lens, const std::size_t c, const F& f, const V& val, const Token& token)
{
    static_assert(!orient_ || is_map_, "Orientations go into the values.");

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

    // The user's count of keys not yet flushed into `approx_sz` is `pending + (ins_c - flushed)`: `pending` from its
    // earlier updates, the rest from this one. It is flushed, and a resize checked for, once `ins_c` reaches `due`.
    std::size_t pending = checkpoint[token.id].unwrap(), flushed = 0;
    const auto step = checkpoint_step.load(std::memory_order_relaxed);
    std::size_t due = (pending < step ? step - pending : 1);

    // Flushes the user's count. Returns whether the table is due for a resize.
    const auto flush_count = [&]() -> bool
    {
        const auto sz = (approx_sz += pending + (ins_c - flushed));
        pending = 0, flushed = ins_c, due = ins_c + checkpoint_step.load(std::memory_order_relaxed);
        return sz >= resize_th.load(std::memory_order_relaxed);
    };

    // Resizes the table, iff it is due.
    const auto resize_now = [&]()
    {
        release();
        if constexpr(mt_)   table_lock.unlock_shared(token.id);
        try_resize();
        if constexpr(mt_)   table_lock.lock_shared(token.id);
        in_run = false;
    };

    // Grows the overflow table, as `try_grow_ov(ov, full)` does.
    const auto grow_ov_now = [&](const bool full)
    {
        const auto* const ov_seen = ov;
        release();
        if constexpr(mt_)   table_lock.unlock_shared(token.id);
        try_grow_ov(ov_seen, full);
        if constexpr(mt_)   table_lock.lock_shared(token.id);
        in_run = false;
    };

    // Accounts for a k-mer inserted, resizing the table if it is due.
    const auto inserted = [&]()
    {
        if(KACHE_HASH_UNLIKELY(++ins_c == due) && flush_count())
            resize_now();
    };

    // An update's lambdas, as the canonical ones but with the k-mer's orientation `fwd` iff `orient_`, which the
    // values are computed from: each branch has its own, as the canonical ones call `f(v)`.
    if constexpr(orient_)
    {
        const auto place_o = [&](const Kmer_Desc& d, const bool fwd, const std::size_t i, const std::size_t j)
        {
            place_at(flat_t(d.key, f(val, fwd)), d.c, i, j);
        };

        const auto found_o = [&](const bool fwd, const std::size_t i, const std::size_t j)
        {
            T.set_val(i, j, f(T.val(i, j), fwd));
        };

        const auto apply_overflow_o = [&](const Kmer_Desc& d, const bool fwd) -> bool
        {
            bool full;
            const auto absent = ov_update(d.key, [&](const mapped_t v){ return f(v, fwd); }, f(val, fwd), full);
            if(KACHE_HASH_UNLIKELY(full))
                return grow_ov_now(true), false;

            if(absent)
            {
                inserted();
                if(add_ov_key(token.id))
                    grow_ov_now(false);
            }

            return true;
        };

        const auto apply_o = [&](const Kmer_Desc& d, const bool fwd) -> bool
        {
            if(!in_run || d.h != run_h)
            {
                release();
                if constexpr(mt_)
                    if(KACHE_HASH_UNLIKELY(table_lock.writer_pending()))  // Yield to a resize, holding no bucket locks.
                        table_lock.unlock_shared(token.id), table_lock.lock_shared(token.id);

                in_run = true, run_h = d.h, spilled = false;
                b = bucket_of(d.h);
            }

            if(!spilled)
            {
                if(held == 0)
                    lock(b), held = 1, cs_b = load_cs(md(b));

                uint32_t e;
                if(const auto j = probe(d.key, d.c, b, cs_b, e); j < B)
                    return found_o(fwd, b, j), true;

                if(e < B)
                {
                    place_o(d, fwd, b, e);
                    cs_b = set_cs(cs_b, e, d.c);
                    T.prefetch_frontier(b, e);
                    return inserted(), true;
                }

                mark_spill(md(b), d.h);
                unlock(b), held = 0;
                std::tie(p, q) = secondaries(d.h);
                spilled = true;
            }

            if(held == 0)
                lock_sorted(b, p, q), held = 3;

            uint32_t e_b, e_p, e_q;
            if(const auto j = probe(d.key, d.c, b, e_b); j < B)
                return found_o(fwd, b, j), true;
            if(const auto j = probe(d.key, d.c, p, e_p); j < B)
                return found_o(fwd, p, j), true;
            if(const auto j = probe(d.key, d.c, q, e_q); j < B)
                return found_o(fwd, q, j), true;

            if(e_p < e_q)
                return place_o(d, fwd, p, e_p), inserted(), true;
            if(e_q < B)
                return place_o(d, fwd, q, e_q), inserted(), true;

            release();
            return apply_overflow_o(d, fwd);
        };

        if constexpr(mt_)   table_lock.lock_shared(token.id);

        // A call's stretches go through one pipeline of the batch describer; once the slot prefetches cover their
        // insertions, it gains build performance. In the lanes, a lane's last minimizer run in a chunk continues into the
        // next chunk's k-mers of the lane, so that chunk boundaries don't split runs, and their locks.
        Sec_Cache sc_slot;
        stream<use_lanes, true, false>(seqs, lens, c, cap_,
            [&](const Kmer_Desc& d){ prefetch_meta<true>(d); },
            [&](const auto& make){ prefetch_slot<true>(make(), sc_slot); },
            [&](const Kmer_Desc& d, std::size_t, std::size_t, const bool fwd)
            {
                while(!apply_o(d, fwd))
                    ;
            });
    }
    else
    {
        // Places the k-mer `d` into the `i`'th bucket at slot `j`.
        const auto place = [&](const Kmer_Desc& d, const std::size_t i, const std::size_t j)
        {
            if constexpr(is_set_)
                place_at(d.key, d.c, i, j);
            else
                place_at(flat_t(d.key, val), d.c, i, j);
        };

        // Applies the update to the value of the present k-mer at the `i`'th bucket's slot `j`.
        const auto found = [&](const std::size_t i, const std::size_t j)
        {
            if constexpr(is_map_)
                T.set_val(i, j, f(T.val(i, j)));
            else
                (void)i, (void)j;
        };

        // Applies the overflow table's update with the k-mer `d`, holding no lock. Returns `false` iff it found the
        // overflow table full, having grown it instead; the update is to be retried then.
        const auto apply_overflow = [&](const Kmer_Desc& d) -> bool
        {
            bool full;
            const auto absent = ov_update(d.key, f, val, full);
            if(KACHE_HASH_UNLIKELY(full))
                return grow_ov_now(true), false;

            if(absent)
            {
                inserted();
                if(add_ov_key(token.id))
                    grow_ov_now(false);
            }

            return true;
        };

        // Applies the update with the k-mer `d`. Returns `false` iff it found the overflow table full, having grown it
        // instead; the update is to be retried then.
        const auto apply = [&](const Kmer_Desc& d) -> bool
        {
            if(!in_run || d.h != run_h)
            {
                release();
                if constexpr(mt_)
                    if(KACHE_HASH_UNLIKELY(table_lock.writer_pending()))  // Yield to a resize, holding no bucket locks.
                        table_lock.unlock_shared(token.id), table_lock.lock_shared(token.id);

                in_run = true, run_h = d.h, spilled = false;
                b = bucket_of(d.h);
            }

            if(!spilled)
            {
                if(held == 0)
                    lock(b), held = 1, cs_b = load_cs(md(b));

                uint32_t e;
                if(const auto j = probe(d.key, d.c, b, cs_b, e); j < B)
                    return found(b, j), true;

                if(e < B)
                {
                    place(d, b, e);
                    cs_b = set_cs(cs_b, e, d.c);

                    // A run's insertions fill its bucket's lines in order, and the slot prefetches predict the run's
                    // first only: placing a key at a line's first slot prefetches the next line.
                    T.prefetch_frontier(b, e);

                    return inserted(), true;
                }

                // The spill is marked before placing beyond, so a lookup finding no mark precedes the placements.
                mark_spill(md(b), d.h);
                unlock(b), held = 0;
                std::tie(p, q) = secondaries(d.h);
                spilled = true;
            }

            if(held == 0)
                lock_sorted(b, p, q), held = 3;

            uint32_t e_b, e_p, e_q;
            if(const auto j = probe(d.key, d.c, b, e_b); j < B)
                return found(b, j), true;
            if(const auto j = probe(d.key, d.c, p, e_p); j < B)
                return found(p, j), true;
            if(const auto j = probe(d.key, d.c, q, e_q); j < B)
                return found(q, j), true;

            if(e_p < e_q)
                return place(d, p, e_p), inserted(), true;
            if(e_q < B)
                return place(d, q, e_q), inserted(), true;

            release();
            return apply_overflow(d);
        };

        if constexpr(mt_)   table_lock.lock_shared(token.id);

        // A call's stretches go through one pipeline of the batch describer; once the slot prefetches cover their
        // insertions, it gains build performance. In the lanes, a lane's last minimizer run in a chunk continues into the
        // next chunk's k-mers of the lane, so that chunk boundaries don't split runs, and their locks.
        Sec_Cache sc_slot;
        stream<use_lanes, true, false>(seqs, lens, c, cap_,
            [&](const Kmer_Desc& d){ prefetch_meta<true>(d); },
            [&](const auto& make){ prefetch_slot<true>(make(), sc_slot); },
            [&](const Kmer_Desc& d, std::size_t, std::size_t)
            {
                while(!apply(d))
                    ;
            });

    }

    release();
    if constexpr(mt_)   table_lock.unlock_shared(token.id);
    checkpoint[token.id].unwrap() = pending + (ins_c - flushed);

    return ins_c;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename G>
__attribute__((flatten)) inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::find(const char* const* const seqs, const std::size_t* const lens, const std::size_t c, const G& g) const
{
    // The current minimizer run, its primary bucket `b`, and its secondary buckets `p` and `q` once computed.
    uint64_t run_h = 0;
    bool sec = false;
    std::size_t b = 0, p = 0, q = 0;

    // Returns what `find` returns for the present k-mer at the `i`'th bucket's slot `j`.
    const auto present = [&](const std::size_t i, const std::size_t j) -> find_ret_t
    {
        (void)i, (void)j;
        if constexpr(is_map_)   return T.val(i, j);
        else                    return true;
    };

    // Whether `g` takes the k-mers' orientations; one that takes either form gets none.
    constexpr bool orient_ = !std::is_invocable_v<const G&, std::size_t, std::size_t, const find_ret_t&> &&
                             std::is_invocable_v<const G&, std::size_t, std::size_t, const find_ret_t&, bool>;

    // Positive queries gain from prefetching their slots, and negative ones lose the work, so the slots are prefetched
    // while the last query hit; queries tend to come in streaks. A call's stretches go through one pipeline, which a
    // pipeline per stretch, such as per read, would restart. Lane runs don't continue across chunks for queries, which
    // take no locks: that costs them.
    bool hit = true;
    Sec_Cache sc_slot;
    if constexpr(!orient_)
        stream<use_lanes, false, true>(seqs, lens, c, cap_,
            [&](const Kmer_Desc& d){ prefetch_meta<false>(d); },
            Lookup_Pre_Slot<false>{*this, sc_slot, hit},
            [&](const Kmer_Desc& d, const std::size_t seq_i, const std::size_t i)
            {
                // NB: Branch-free, as a branch on run starts mispredicts once per run.
                b = bucket_of(d);
                sec &= (d.h == run_h);
                run_h = d.h;

                // Absent past a full primary bucket too, unless its spill filter marks the minimizer.
                bool full;
                if(const auto j = probe(d.key, d.c, b, full); j < B)
                    return hit = true, g(seq_i, i, present(b, j));
                if(!full || !may_spill(md(b), d.h))
                    return hit = false, g(seq_i, i, find_ret_t(null_val));

                if(!sec)
                    std::tie(p, q) = secondaries(d.h), sec = true;

                bool full_p, full_q;
                if(const auto j = probe(d.key, d.c, p, full_p); j < B)
                    return hit = true, g(seq_i, i, present(p, j));
                if(const auto j = probe(d.key, d.c, q, full_q); j < B)
                    return hit = true, g(seq_i, i, present(q, j));

                const auto r = (full_p && full_q) ? find_in_overflow(d.key) : find_ret_t(null_val);
                hit = static_cast<bool>(r);
                g(seq_i, i, r);
            });
    else    // As above, with the orientations.
        stream<use_lanes, false, true>(seqs, lens, c, cap_,
            [&](const Kmer_Desc& d){ prefetch_meta<false>(d); },
            Lookup_Pre_Slot<false>{*this, sc_slot, hit},
            [&](const Kmer_Desc& d, const std::size_t seq_i, const std::size_t i, const bool fwd)
            {
                b = bucket_of(d);
                sec &= (d.h == run_h);
                run_h = d.h;

                bool full;
                if(const auto j = probe(d.key, d.c, b, full); j < B)
                    return hit = true, g(seq_i, i, present(b, j), fwd);
                if(!full || !may_spill(md(b), d.h))
                    return hit = false, g(seq_i, i, find_ret_t(null_val), fwd);

                if(!sec)
                    std::tie(p, q) = secondaries(d.h), sec = true;

                bool full_p, full_q;
                if(const auto j = probe(d.key, d.c, p, full_p); j < B)
                    return hit = true, g(seq_i, i, present(p, j), fwd);
                if(const auto j = probe(d.key, d.c, q, full_q); j < B)
                    return hit = true, g(seq_i, i, present(q, j), fwd);

                const auto r = (full_p && full_q) ? find_in_overflow(d.key) : find_ret_t(null_val);
                hit = static_cast<bool>(r);
                g(seq_i, i, r, fwd);
            });
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::find(const char* const kmer) const -> find_ret_t
{
    if(DNA_Utility::acgt_prefix(kmer, k) < k)
        return null_val;

    // The k-mer's descriptor, as the describers give it.
    const auto key = Kmer<k>(kmer).canonical();
    const auto x = key.as_int();
    uint64_t h;
    Kmer_Lanes<k, l>::minimizer_hashes(&x, 1, &h);
    const auto c = kmer_checksum(x);

    // Returns what `find` returns for the present k-mer at the `i`'th bucket's slot `j`.
    const auto present = [&](const std::size_t i, const std::size_t j) -> find_ret_t
    {
        (void)i, (void)j;
        if constexpr(is_map_)   return T.val(i, j);
        else                    return true;
    };

    // Probed as the sequence `find`s probe: absent past a full primary bucket too, unless its spill filter marks the
    // minimizer.
    const auto b = bucket_of(h);
    bool full;
    if(const auto j = probe(key, c, b, full); j < B)
        return present(b, j);
    if(!full || !may_spill(md(b), h))
        return null_val;

    const auto [p, q] = secondaries(h);
    bool full_p, full_q;
    if(const auto j = probe(key, c, p, full_p); j < B)
        return present(p, j);
    if(const auto j = probe(key, c, q, full_q); j < B)
        return present(q, j);

    return (full_p && full_q) ? find_in_overflow(key) : find_ret_t(null_val);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename F>
__attribute__((flatten)) inline std::size_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::update(const char* const* const seqs, const std::size_t* const lens, const std::size_t c, const F& f) requires (is_map_)
{
    // The current minimizer run, its primary bucket `b`, and its secondary buckets `p` and `q` once computed.
    uint64_t run_h = 0;
    bool sec = false;
    std::size_t b = 0, p = 0, q = 0;

    std::size_t found = 0;  // Number of k-mers present.

    // The k-mers are looked up as in `find`, with the slots prefetched for writing.
    bool hit = true;
    Sec_Cache sc_slot;
    stream<use_lanes, false, true>(seqs, lens, c, cap_,
        [&](const Kmer_Desc& d){ prefetch_meta<false>(d); },
        Lookup_Pre_Slot<true>{*this, sc_slot, hit},
        [&](const Kmer_Desc& d, std::size_t, std::size_t)
        {
            b = bucket_of(d);
            sec &= (d.h == run_h);
            run_h = d.h;

            bool full;
            if(const auto j = probe(d.key, d.c, b, full); j < B)
                return hit = true, found++, T.update_val(b, j, f);
            if(!full || !may_spill(md(b), d.h))
                return (void)(hit = false);

            if(!sec)
                std::tie(p, q) = secondaries(d.h), sec = true;

            bool full_p, full_q;
            if(const auto j = probe(d.key, d.c, p, full_p); j < B)
                return hit = true, found++, T.update_val(p, j, f);
            if(const auto j = probe(d.key, d.c, q, full_q); j < B)
                return hit = true, found++, T.update_val(q, j, f);

            hit = (full_p && full_q && update_in_overflow(d.key, f));
            found += hit;
        });

    return found;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline std::size_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::count(const char* const* const seqs, const std::size_t* const lens, const std::size_t c) const
{
    std::size_t hit_c = 0;
    find(seqs, lens, c, [&](std::size_t, std::size_t, const find_ret_t& r){ hit_c += static_cast<bool>(r); });
    return hit_c;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename F>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::for_each(const F& f, const std::size_t w_c, const std::size_t w_id) const
{
    assert(w_id < w_c);

    // The worker's range of buckets, then its range of the overflow table.
    for(std::size_t b = cap_ * w_id / w_c; b < cap_ * (w_id + 1) / w_c; ++b)
        for(std::size_t j = 0, n = bucket_size(b); j < n; ++j)
            if constexpr(is_set_)
                f(T.key(b, j));
            else
                f(T.key(b, j), T.val(b, j));

    auto it = ov->iterator(w_c, w_id);
    Kmer<k> key;
    ov_val_t val;
    while(it.next(key, val))
        if constexpr(is_set_)
            f(key);
        else
            f(key, val);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::key_distribution() const -> std::tuple<double, double, double>
{
    std::atomic<std::size_t> p = 0, s = 0;
    std::vector<std::thread> W;
    W.reserve(resize_worker_c);
    for(uint64_t w = 0; w < resize_worker_c; ++w)
        W.emplace_back([&](const uint64_t w_id)
        {
            std::size_t p_w = 0, s_w = 0;
            for(std::size_t b = cap_ * w_id / resize_worker_c; b < cap_ * (w_id + 1) / resize_worker_c; ++b)
            {
                const auto d = distribution(b);
                p_w += d.first, s_w += d.second;
            }

            p += p_w, s += s_w;
        }, w);

    std::for_each(W.begin(), W.end(), [](auto& w){ w.join(); });

    const auto o = overflow_size();
    const auto sz = p + s + o;
    return {p * 100.0 / sz, s * 100.0 / sz, o * 100.0 / sz};
}


}



#endif
