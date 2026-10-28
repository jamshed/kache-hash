
#ifndef STREAMING_KMER_HASHSET_HPP
#define STREAMING_KMER_HASHSET_HPP



#include "DNA_Utility.hpp"
#include "Kmer.hpp"
#include "Kmer_Utility.hpp"
#include "Kmer_Batch.hpp"
#include "Kmer_Lanes.hpp"
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
#include <iostream>
#include <vector>
#include <cstdlib>
#include <thread>
#include <algorithm>
#include <barrier>
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
    static constexpr double of_default = 1.0 / 32;  // Initial overflow factor: the overflow table starts with room for
                                                    // this fraction of the main table's capacity, and grows on its own.
    static constexpr double ov_lf = 0.75;   // Load factor at which the overflow table grows.
    static constexpr std::size_t ov_min_cap = 1024; // Fewest keys the overflow table has room for.
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

    // Metadata of the k-mers in a bucket. The top bits of `min_coord[1..B)`
    // are its spill filter: they mark the minimizers whose keys left it, their
    // primary bucket.
    struct Metadata
    {
        uint8_t cs[B];  // Checksums of the k-mers.
        uint8_t min_coord[B];   // Minimizer coordinates within the k-mers.

        static constexpr uint8_t coord_mask = 0b0111'1111;  // Coordinate bits of a `min_coord`.
        static constexpr uint8_t lock_mask = 0b1000'0000;   // Bit of `min_coord[0]` used to lock the bucket.
        static constexpr uint8_t spill_mask = 0b1000'0000;  // Bit of `min_coord[1..B)` used by the spill filter.
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


    // Storage of a table's buckets: each bucket holds the metadata of its
    // keys, then their slots, so that an operation finds both on one page. The
    // buckets sit on huge pages where supported.
    class Store
    {
        // A bucket: the metadata of its keys, then their slots.
        struct alignas(L1_CACHE_LINE_SIZE) Bucket
        {
            Metadata md;    // Metadata of the keys.
            flat_t slot[B]; // The keys, with their values in a map.
        };

        Bucket* T = nullptr;    // The buckets.

    public:

        // Returns the metadata of the `i`'th bucket.
        Metadata& md(const std::size_t i) const { return T[i].md; }

        // Returns the `j`'th slot of the `i`'th bucket.
        flat_t& slot(const std::size_t i, const std::size_t j) const { return T[i].slot[j]; }

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

    // Returns the `j`'th slot of the `i`'th bucket.
    flat_t& slot(const std::size_t i, const std::size_t j) const { return T.slot(i, j); }

    // Returns the metadata of the `i`'th new bucket in a resize.
    Metadata& md_new(const std::size_t i) const { return T_new.md(i); }

    // Returns the `j`'th slot of the `i`'th new bucket in a resize.
    flat_t& slot_new(const std::size_t i, const std::size_t j) const { return T_new.slot(i, j); }

    std::atomic_uint16_t registered_thread_count;   // Number of threads registered to use the table.

    std::atomic_uint64_t approx_sz; // Size of the hash table, short of the user-local counts not yet flushed into it.
    std::vector<Padded<uint64_t>> checkpoint;   // `checkpoint[i]` is the number of keys the `i`'th user has added and not yet flushed into `approx_sz`.
    static constexpr uint64_t resize_checkpoint = 16384;    // Most keys a user adds before flushing its count and checking for a resize.
    std::atomic_uint64_t checkpoint_step;   // Number of keys a user adds before flushing its count and checking for a resize: unflushed keys stay under an eighth of the resize threshold.
    RW_Lock<max_user> table_lock; // Readers-writers lock to resize the table exclusively.

    const uint64_t resize_worker_c; // Number of worker threads to spawn for a concurrent resize.

    std::atomic_uint64_t ov_sz; // Size of the overflow table, short of the user-local counts not yet flushed into it.
    std::vector<Padded<uint64_t>> ov_checkpoint;    // `ov_checkpoint[i]` is the number of keys the `i`'th user has added to the overflow table and not yet flushed into `ov_sz`.
    std::atomic_uint64_t ov_th; // Size of the overflow table at which it grows.
    std::atomic_uint64_t ov_step;   // Number of keys a user adds to the overflow table before flushing its count and checking for growth.


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

    // Returns the index in `[1, B)` of the coordinate whose top bit marks the
    // minimizer hash `h` in a spill filter.
    // NB: From bits 32–36 of `h`: the bucket takes its low bits, and window
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

    // Directly places the element `e` with checksum `c` and minimizer-
    // coordinate `m` into the `i`'th bucket at index `j`. The bucket's lock-bit
    // in `min_coord[0]` is kept as is.
    void place_at(const flat_t& e, uint8_t c, uint8_t m, std::size_t i, std::size_t j);

    // Directly places the element `e` with checksum `c` and minimizer-
    // coordinate `m` into the `i`'th bucket at index `j` of the new resized
    // table. The bucket's lock-bit in `min_coord[0]` is kept as is.
    void place_at_resize(const flat_t& e, uint8_t c, uint8_t m, std::size_t i, std::size_t j);

    // Tries to insert the element `x` with checksum `c`, minimizer-coordinate
    // `m`, and minimizer hash `h` into its primary bucket `i` of the newly
    // resized table. Returns whether there was empty space; if not, marks `h`
    // in the bucket's spill filter.
    bool try_insert_at_resize(const flat_t& x, uint8_t c, uint8_t m, uint64_t h, std::size_t i);

    // Inserts the element `x` with checksum `c`, minimizer-coordinate `m`, and
    // minimizer hash `h` into the resized table. Returns whether it goes to the
    // overflow table.
    bool insert_at_resize(const flat_t& x, uint8_t c, uint8_t m, uint64_t h);

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
    // and a key at the old bucket `b` stays there or moves to `b + old_cap`,
    // if it is at its primary bucket; else it is set aside and re-inserted.
    void resize();

    // Descriptor of a k-mer for the streaming operations: what the table
    // operations need of it.
    struct Kmer_Desc
    {
        Kmer<k> key;    // Canonical form of the k-mer.
        uint64_t h; // Hash of its minimizer.
        uint8_t c;  // Its checksum.
        uint8_t m;  // Coordinate of its minimizer, oriented per the canonical form; 0 if not described.
    };

    static constexpr std::size_t line_slots = std::max<std::size_t>(1, L1_CACHE_LINE_SIZE / sizeof(flat_t));  // Number of slots in a cache line.
    static constexpr std::size_t slot_dist = 16;    // How many k-mers ahead of their table operations the streaming
                                                    // operations prefetch the slots of their keys in their primary
                                                    // buckets, or else what they touch beyond. In the lanes, that
                                                    // prefetches all of a chunk's slots a chunk ahead.

    // Applies `op(d, q, i)` to the descriptor `d` of each k-mer of the
    // sequences `seqs[0..c)` of lengths `lens[0..c)`, the k-mer being at index
    // `i` of the `q`'th sequence. Non-ACGT characters split the sequences. A
    // stretch's k-mers come in order, but stretches may interleave. Applies
    // `pre_meta(d)` to each descriptor some k-mers ahead of `op`, and
    // `pre_slot(make)` `slot_dist` k-mers ahead, where `make()` returns the
    // descriptor. Skipping a prefetch then skips making its descriptor. With
    // `lanes_`, the k-mers are described in SIMD lanes, and iff `cont_`, lane
    // runs continue across chunks; else iff `pipe_`, the stretches of ACGT
    // bases go through one pipeline of the batch describer, else a pipeline
    // each. The minimizer coordinates are described iff `coords_`. If
    // `pre_slot` has `active()`, the lanes skip a chunk's slot prefetches while
    // it is `false`.
    // NB: `pre_meta` takes every descriptor, not just each minimizer run's
    // first: a prefetch per k-mer costs less than branching on run starts.
    template <bool lanes_, bool pipe_, bool cont_, bool coords_, typename Pre_Meta_, typename Pre_Slot_, typename Op_>
    static void stream(const char* const* seqs, const std::size_t* lens, std::size_t c, const Pre_Meta_& pre_meta, const Pre_Slot_& pre_slot, const Op_& op);

    // Applies `stream`'s callbacks to the k-mers of the sequences `seqs[0..c)`
    // of lengths `lens[0..c)` through the SIMD lanes: long stretches split over
    // the lanes on their own, short ones gathered a lane each.
    template <bool cont_, bool coords_, typename Pre_Meta_, typename Pre_Slot_, typename Op_>
    static void stream_lanes_all(const char* const* seqs, const std::size_t* lens, std::size_t c, const Pre_Meta_& pre_meta, const Pre_Slot_& pre_slot, const Op_& op);

    // Applies `stream`'s callbacks to the `n` k-mers of the ACGT stretch `t`,
    // whose first k-mer has index `base` in the `q`'th sequence, describing
    // them in chunks: the metadata of a chunk's buckets is prefetched one to
    // two chunks ahead.
    template <bool coords_, typename Pre_Meta_, typename Pre_Slot_, typename Op_>
    static void stream_batched(const char* t, std::size_t n, std::size_t q, std::size_t base, const Pre_Meta_& pre_meta, const Pre_Slot_& pre_slot, const Op_& op);

    static constexpr std::size_t lane_c = Kmer_Lanes<k, l>::L;  // Number of SIMD lanes describing stretches at once.
    static constexpr std::size_t lane_min_kmers = 4096; // Fewest k-mers in a stretch for it to be split over the lanes on its own.
    static constexpr std::size_t lane_depth = 2;    // Number of lane chunks in flight: the one operated on, and the one
                                                    // ahead, whose metadata is prefetched once described.
    static constexpr std::size_t lane_group_c = 16; // Most groups of short stretches, such as reads, that go through one
                                                    // pipeline of lane chunks: 1 restarts the pipeline per group.

    // The streaming operations describe k-mers in SIMD lanes only for targets
    // with AVX-512, whose cores run 256-bit vectors at full width. Cores that
    // split them into 128-bit halves, such as Sierra Forest's E-cores, run the
    // batch describer faster.
#ifdef __AVX512F__
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
    template <bool coords_, typename Pre_Meta_, typename Pre_Slot_, typename Op_>
    static void stream_stretches(const Stretch* S, std::size_t c, const Pre_Meta_& pre_meta, const Pre_Slot_& pre_slot, const Op_& op);

    // Applies `stream`'s callbacks to the k-mers of `g_c` groups of
    // stretches: the `g`'th group's are `S[g * lane_c + j]` for `j` in
    // `[0, G_c[g])`, with `G_c[g] <= lane_c`. Each group goes through the SIMD
    // lanes at once, a stretch per lane, and the groups through one pipeline of
    // chunks. The metadata of a chunk's buckets is prefetched once it is
    // described, before the operations on the chunk ahead of it. A chunk's
    // k-mers go lane by lane, a lane's in order. All of a chunk's slots are
    // prefetched a chunk ahead. Iff `cont_`, lane runs continue across chunks.
    template <bool cont_, bool coords_, typename Pre_Meta_, typename Pre_Slot_, typename Op_>
    static void stream_lanes(const Stretch* S, const std::size_t* G_c, std::size_t g_c, const Pre_Meta_& pre_meta, const Pre_Slot_& pre_slot, const Op_& op);

    // Describes the next chunk of the lanes `kl` into `out` without minimizer
    // coordinates, as lookups do. Returns what `Kmer_Lanes::next` does.
    // NB: Out of line, so that the lanes' unrolled body stays apart from the
    // lookups' loops. The batch describer runs faster inlined.
    __attribute__((noinline)) static bool describe_lookups(Kmer_Lanes<k, l>& kl, typename Kmer_Lanes<k, l>::Chunk& out);

    // Returns the index of the slot holding the k-mer `key` with checksum `c`
    // and minimizer coordinate `m` in the `i`'th bucket, or `B` if it is absent
    // there. Puts the index of the bucket's first empty slot at `e`, `B` if it
    // is full.
    // NB: `m` goes unused. Passing it anyway steers GCC's register allocation
    // in the flattened set queries away from a spill per k-mer.
    uint32_t probe(Kmer<k> key, uint8_t c, uint8_t m, std::size_t i, uint32_t& e) const { return probe(key, c, m, i, load_cs(md(i)), e); }

    // Returns the index of the slot holding the k-mer `key` with checksum `c`
    // and minimizer coordinate `m` in the `i`'th bucket, which has the
    // checksums `c_vec`; or `B` if it is absent there. Puts the index of the
    // bucket's first empty slot at `e`, `B` if it is full.
    uint32_t probe(Kmer<k> key, uint8_t c, uint8_t m, std::size_t i, cs_vec c_vec, uint32_t& e) const;

    // Prefetches the metadata of the bucket of the k-mer `d`, with write
    // intent iff `write_`.
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
    // iff `write_`, with write intent iff `write_`. If the operation goes
    // beyond the bucket, prefetches its secondary buckets' metadata instead,
    // computing them through the cache `sc`, and its overflow-table slot. A
    // lookup goes beyond only if the spill filter marks its minimizer.
    template <bool write_> void prefetch_slot(const Kmer_Desc& d, Sec_Cache& sc) const;

    // Prefetches what an operation on the k-mer `key` with minimizer hash `h`
    // touches beyond its full primary bucket, as `prefetch_slot` does,
    // computing the secondary buckets through the cache `sc`.
    // NB: Out of line, as inlined it bloats the operations' loops by 12-18
    // instructions per k-mer. It takes the k-mer and the hash by value, as a
    // descriptor by reference would be kept in memory for the call.
    template <bool write_> __attribute__((noinline, cold)) void prefetch_spill(Kmer<k> key, uint64_t h, Sec_Cache& sc) const;

    // Lookups' slot prefetches, as `stream` takes them: on while the last
    // lookup hit. The lanes ask once per chunk whether they are on.
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

    // Inserts the k-mer `key` with checksum `c` and coordinate `m` into the
    // overflow table iff it is absent there, with value `val` in a map; also in
    // a map, updates the value `v` of a present one to `f(v)`. Returns whether
    // it was absent. Sets `full` iff it is absent and finds no room; then
    // nothing changes.
    // NB: Out of line, as the streaming updates rarely reach the overflow
    // table.
    template <typename F, typename V>
    __attribute__((noinline)) bool ov_update(Kmer<k> key, uint8_t c, uint8_t m, const F& f, const V& val, bool& full);

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
    set_ov_th();

    // Map the table up front, in parallel, rather than on first touches by the updates. The metadata comes zeroed:
    // checksum `0` denotes absence of a key.
    T.allocate(cap_, resize_worker_c);
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
inline std::size_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::filter_pos(const uint64_t h)
{
    const auto x = (h >> 32) & (B - 1);
    return x + (x == 0);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::may_spill(const Metadata& md, const uint64_t h)
{
    return md.min_coord[filter_pos(h)] & Metadata::spill_mask;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::mark_spill(Metadata& md, const uint64_t h)
{
    md.min_coord[filter_pos(h)] |= Metadata::spill_mask;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::clear_spills(Metadata& md)
{
    for(std::size_t j = 1; j < B; ++j)
        md.min_coord[j] &= ~Metadata::spill_mask;
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
            const auto h = minimizer_hash(key(slot(b, j)), md(b).min_coord[j] & Metadata::coord_mask);
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
                // A waiter spins on loads, pausing, and tries the CAS only once the bucket looks unlocked: spinning
                // CASes would take the metadata line from the holder, which works on it; builds collapsed on hot
                // buckets at 56–128 threads.
                auto& w = (resize_ ? md_new(b) : md(b)).min_coord[0];
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
inline bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::try_insert_at_resize(const flat_t& e, const uint8_t c, const uint8_t m, const uint64_t h, const std::size_t i)
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

    place_at_resize(e, c, m, i, j);
    unlock<true>(i);
    return true;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::insert_at_resize(const flat_t& x, const uint8_t c, const uint8_t m, const uint64_t h)
{
    const auto idx_mask = cap_ - 1;

    const auto b = h & idx_mask;
    if(try_insert_at_resize(x, c, m, h, b))
        return false;

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

        bool full;
        if constexpr(is_set_)
            ov_new->insert(x, Overflow_Set_Val(c, m), full);
        else
            ov_new->insert(x.first, Overflow_Map_Val(c, m, x.second), full);

        if(full)
        {
            std::cerr << "Overflow table of capacity " << ov_new->capacity() << " is full in a resize. Aborting.\n";
            std::abort();
        }

        return true;
    }

    return false;
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
                typedef std::conditional_t<!is_map_, Overflow_Set_Val, Overflow_Map_Val> ov_val_t;
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
    if(!T.grow(old_cap, 2 * old_cap))   // The new buckets come zeroed: checksum `0` denotes absence of a key.
    {
        std::cerr << "Growing the table to " << 2 * old_cap * B << " slots failed. Aborting.\n";
        std::abort();
    }

    cap_ *= 2;
    T_new = T;  // The resize-insertion helpers act on this.
    ov_new = new ov_t(new_ov_size(), resize_worker_c);
    std::atomic<uint64_t> ov_new_c{0};  // Number of keys moved into the new overflow table.

    // In the first phase, a worker owns its range of old buckets `b` and their new halves `b + old_cap`: the keys
    // at their primary buckets stay, compacted, or move to the new half, without locks. The keys at secondary
    // buckets are set aside, and with the overflow table's, re-inserted with locks in the second phase.
    const auto worker_c = resize_worker_c;
    std::barrier phase_end(static_cast<std::ptrdiff_t>(worker_c));
    std::vector<std::thread> W;
    W.reserve(worker_c);
    for(uint64_t i = 0; i < worker_c; ++i)
        W.emplace_back([&](const auto w_id)
        {
            struct Held { flat_t x; uint8_t c, m; };
            std::vector<Held> held; // This worker's keys at secondary buckets.

            const auto b_s = old_cap * w_id / worker_c, b_e = old_cap * (w_id + 1) / worker_c;
            for(std::size_t b = b_s; b < b_e; ++b)
            {
                const auto nb = b + old_cap;
                std::size_t lo = 0, hi = 0; // Numbers of keys kept at `b`, and moved to `nb`.
                const auto n = bucket_size(b);
                for(std::size_t j = 0; j < n; ++j)
                {
                    const uint8_t c = md(b).cs[j];
                    const uint8_t m = md(b).min_coord[j] & Metadata::coord_mask;
                    const auto h = minimizer_hash(key(slot(b, j)), m);
                    if((h & (old_cap - 1)) != b)
                    {
                        held.push_back({slot(b, j), c, m});
                        continue;
                    }

                    if(h & old_cap)
                        slot(nb, hi) = slot(b, j), md(nb).cs[hi] = c, md(nb).min_coord[hi] = m, hi++;
                    else
                    {
                        if(lo != j)
                            slot(b, lo) = slot(b, j), md(b).cs[lo] = c, md(b).min_coord[lo] = m;
                        lo++;
                    }
                }

                for(std::size_t j = lo; j < n; ++j)
                    md(b).cs[j] = 0, md(b).min_coord[j] = 0;

                clear_spills(md(b));    // The second phase marks the spills anew.
            }

            phase_end.arrive_and_wait();

            // The rest go through a pipeline: a key's new primary bucket is prefetched `pipe_sz` keys ahead of its
            // insertion.
            constexpr std::size_t pipe_sz = 16;
            struct Pending { flat_t x; uint64_t h; uint8_t c, m; };
            Pending pipe[pipe_sz];
            std::size_t pipe_c = 0;
            uint64_t ov_n = 0;  // Number of keys this worker moves into the new overflow table.
            const auto push = [&](const flat_t& x, const uint8_t c, const uint8_t m)
            {
                const auto h = minimizer_hash(key(x), m);
                const auto& meta = md(h & (cap_ - 1));
                __builtin_prefetch(meta.cs, 1, 3);
                if constexpr(sizeof(Metadata) > L1_CACHE_LINE_SIZE)
                    __builtin_prefetch(meta.min_coord, 1, 3);

                auto& slot = pipe[pipe_c % pipe_sz];
                if(pipe_c >= pipe_sz)
                    ov_n += insert_at_resize(slot.x, slot.c, slot.m, slot.h);
                slot = Pending{x, h, c, m};
                pipe_c++;
            };

            for(const auto& e : held)
                push(e.x, e.c, e.m);
            std::vector<Held>().swap(held);

            auto it = ov->iterator(worker_c, w_id);
            typedef std::conditional_t<!is_map_, Overflow_Set_Val, Overflow_Map_Val> ov_val_t;
            Kmer<k> key;
            ov_val_t val;
            while(it.next(key, val))
            {
                if constexpr(is_set_)
                    push(key, val.cs, val.min_coord);
                else
                    push(flat_t(key, val.val), val.cs, val.min_coord);
            }

            for(std::size_t e = (pipe_c > pipe_sz ? pipe_c - pipe_sz : 0); e < pipe_c; ++e)
            {
                const auto& slot = pipe[e % pipe_sz];
                ov_n += insert_at_resize(slot.x, slot.c, slot.m, slot.h);
            }

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
template <bool lanes_, bool pipe_, bool cont_, bool coords_, typename Pre_Meta_, typename Pre_Slot_, typename Op_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::stream(const char* const* const seqs, const std::size_t* const lens, const std::size_t c, const Pre_Meta_& pre_meta, const Pre_Slot_& pre_slot, const Op_& op)
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
                    if constexpr(!pipe_)
                        stream_batched<coords_>(seq + s, e - s - (k - 1), q, s, pre_meta, pre_slot, op);
                    else
                    {
                        S[s_c++] = {seq + s, e - s, q, s};
                        if(s_c == group_sz)
                            stream_stretches<coords_>(S, s_c, pre_meta, pre_slot, op), s_c = 0;
                    }
                }

                s = e + 1;
            }
        }

        if(s_c > 0)
            stream_stretches<coords_>(S, s_c, pre_meta, pre_slot, op);
    }
    else
        stream_lanes_all<cont_, coords_>(seqs, lens, c, pre_meta, pre_slot, op);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool cont_, bool coords_, typename Pre_Meta_, typename Pre_Slot_, typename Op_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::stream_lanes_all(const char* const* const seqs, const std::size_t* const lens, const std::size_t c, const Pre_Meta_& pre_meta, const Pre_Slot_& pre_slot, const Op_& op)
{
    // Groups of stretches gathered to go through the lanes in one pipeline: a long stretch split over the lanes on
    // its own, or short ones a lane each.
    Stretch G[lane_group_c * lane_c];
    std::size_t G_c[lane_group_c];  // Number of stretches in each group.
    std::size_t g_c = 0;    // Number of groups gathered, the last one possibly partial.
    const auto flush = [&]()
    {
        stream_lanes<cont_, coords_>(G, G_c, g_c, pre_meta, pre_slot, op);
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
        stream_batched<coords_>(P[j].t, P[j].len - (k - 1), P[j].q, P[j].base, pre_meta, pre_slot, op);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool coords_, typename Pre_Meta_, typename Pre_Slot_, typename Op_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::stream_batched(const char* const t, const std::size_t n, const std::size_t q, const std::size_t base, const Pre_Meta_& pre_meta, const Pre_Slot_& pre_slot, const Op_& op)
{
    typedef Kmer_Batch<k, l> batch_t;
    batch_t kb;
    typename batch_t::Chunk ch[2];  // The chunk being operated on, and the next; each is refilled two chunks ahead once done.

    const auto desc = [](const typename batch_t::Chunk& c, const std::size_t j){ return Kmer_Desc{Kmer<k>(c.key[j]), c.h[j], c.c[j], coords_ ? c.m[j] : uint8_t(0)}; };
    const auto fill = [&](typename batch_t::Chunk& c)
    {
        kb.template next<coords_>(c);
        for(std::size_t j = 0; j < c.n; ++j)
            pre_meta(desc(c, j));
    };

    kb.reset(t, n);
    fill(ch[0]), fill(ch[1]);
    std::size_t first = base;   // Index of the current chunk's first k-mer.
    for(std::size_t cur = 0; ch[cur].n > 0; cur ^= 1)
    {
        const auto& A = ch[cur];
        const auto& A_next = ch[cur ^ 1];
        for(std::size_t j = 0; j < A.n; ++j)
        {
            const auto x = j + slot_dist;
            if(x < A.n)
                pre_slot([&](){ return desc(A, x); });
            else if(x - A.n < A_next.n)
                pre_slot([&](){ return desc(A_next, x - A.n); });

            op(desc(A, j), q, first + j);
        }

        first += A.n;
        fill(ch[cur]);
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool coords_, typename Pre_Meta_, typename Pre_Slot_, typename Op_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::stream_stretches(const Stretch* const S, const std::size_t c, const Pre_Meta_& pre_meta, const Pre_Slot_& pre_slot, const Op_& op)
{
    typedef Kmer_Batch<k, l> batch_t;
    constexpr std::size_t depth = 3;    // Number of chunks in flight: the one operated on, and the next ones, whose metadata is being prefetched. Short stretches make short chunks, so two ahead keep the prefetches early enough.
    batch_t kb;
    typename batch_t::Chunk ch[depth];  // Each is refilled `depth` chunks ahead once done.
    std::size_t ch_s[depth];    // Stretch of each chunk; chunks don't span stretches.
    std::size_t ch_first[depth];    // Index of each chunk's first k-mer in its stretch.

    std::size_t s_next = 0; // Next stretch to describe.
    std::size_t s_cur = 0;  // Stretch being described.
    std::size_t at = 0; // Number of k-mers of the stretch described.
    bool active = false;    // Whether a stretch is being described.

    const auto desc = [](const typename batch_t::Chunk& x, const std::size_t j){ return Kmer_Desc{Kmer<k>(x.key[j]), x.h[j], x.c[j], coords_ ? x.m[j] : uint8_t(0)}; };
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

            if(kb.template next<coords_>(ch[x]) == 0)
            {
                active = false;
                continue;
            }

            ch_s[x] = s_cur, ch_first[x] = at, at += ch[x].n;
            for(std::size_t j = 0; j < ch[x].n; ++j)
                pre_meta(desc(ch[x], j));
            return;
        }
    };

    for(std::size_t x = 0; x < depth; ++x)
        fill(x);
    for(std::size_t cur = 0; ch[cur].n > 0; cur = (cur + 1) % depth)
    {
        const auto& A = ch[cur];
        const auto& A_next = ch[(cur + 1) % depth];
        const auto& st = S[ch_s[cur]];
        for(std::size_t j = 0; j < A.n; ++j)
        {
            const auto x = j + slot_dist;
            if(x < A.n)
                pre_slot([&](){ return desc(A, x); });
            else if(x - A.n < A_next.n)
                pre_slot([&](){ return desc(A_next, x - A.n); });

            op(desc(A, j), st.q, st.base + ch_first[cur] + j);
        }

        fill(cur);
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool cont_, bool coords_, typename Pre_Meta_, typename Pre_Slot_, typename Op_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::stream_lanes(const Stretch* const S, const std::size_t* const G_c, const std::size_t g_c, const Pre_Meta_& pre_meta, const Pre_Slot_& pre_slot, const Op_& op)
{
    typedef Kmer_Lanes<k, l> lanes_t;
    constexpr auto T_c = lanes_t::T;
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
    const auto desc = [&](const std::size_t x, const std::size_t j, const std::size_t t)
    {
        return Kmer_Desc{Kmer<k>(ch[x].key[t][j]), ch[x].h[t][j], static_cast<uint8_t>(ch[x].c[t][j]), coords_ ? static_cast<uint8_t>(ch[x].m[t][j]) : uint8_t(0)};
    };
    const auto describe = [&](const std::size_t x)
    {
        if constexpr(coords_)   return kl.template next<true>(ch[x]);
        else                    return describe_lookups(kl, ch[x]);
    };
    const auto fill_and_prefetch = [&](const std::size_t x) -> bool
    {
        while(!describe(x)) // The group is done; describe the next one, if any.
        {
            if(g_d + 1 == g_c)
                return false;
            start(++g_d);
        }

        ch_g[x] = g_d;
        const auto c = G_c[g_d];
        if(ch[x].any)
            for(std::size_t j = 0; j < c; ++j)
            {
                const auto [t_b, t_e] = ch[x].range(j);
                for(std::size_t t = t_b; t < t_e; ++t)
                    pre_meta(desc(x, j, t));
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
                            op(desc(cur_next, j, u), S_g[j].q, base_next + u);

                        done[j] = u - u_b;
                    }
            }

        live[cur] = fill_and_prefetch(cur);
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::describe_lookups(Kmer_Lanes<k, l>& kl, typename Kmer_Lanes<k, l>::Chunk& out)
{
    return kl.template next<false>(out);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline uint32_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::probe(const Kmer<k> key, const uint8_t c, const uint8_t m, const std::size_t i, const cs_vec c_vec, uint32_t& e) const
{
    (void)m;
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
template <bool write_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::prefetch_meta(const Kmer_Desc& d) const
{
    const auto& meta = md(d.h & (cap_ - 1));
    __builtin_prefetch(meta.cs, write_, 3);
    if constexpr(write_ && sizeof(Metadata) > L1_CACHE_LINE_SIZE)  // Updates lock and write the coordinates too.
        __builtin_prefetch(meta.min_coord, write_, 3);
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
template <bool write_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::prefetch_slot(const Kmer_Desc& d, Sec_Cache& sc) const
{
    const auto i = d.h & (cap_ - 1);
    const auto c_vec = load_cs(md(i));
    const auto cand = eq_mask(c_vec, splat_cs(d.c));
    if(!cand)
    {
        const auto empty = eq_mask(c_vec, zero_cs());
        if(!empty)  // The operation goes beyond the primary bucket.
        {
            if constexpr(!write_)
                if(!may_spill(md(i), d.h))
                    return;

            return prefetch_spill<write_>(d.key, d.h, sc);
        }

        // Updates predicted to insert prefetch the line after their predicted slot's too: a run's insertions cross
        // into it, and with the insertion-frontier prefetches alone, their stores still filled the store buffer.
        if constexpr(write_)
        {
            const auto e = std::countr_zero(empty);
            __builtin_prefetch(&slot(i, e), write_, 3);
            if(e + line_slots < B)
                __builtin_prefetch(&slot(i, e + line_slots), write_, 3);
        }

        return;
    }

    __builtin_prefetch(&slot(i, std::countr_zero(cand)), write_, 3);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Lookup_Pre_Slot::active() const
{
    return hit;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename Make_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Lookup_Pre_Slot::operator()(const Make_& make) const
{
    if(hit)
        M.template prefetch_slot<false>(make(), sc);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename F, typename V>
bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::ov_update(const Kmer<k> key, const uint8_t c, const uint8_t m, const F& f, const V& val, bool& full)
{
    if constexpr(is_set_)
        return (void)f, (void)val, !ov->insert(key, Overflow_Set_Val(c, m), full);
    else
        return !ov->upsert(key, Overflow_Map_Val(c, m, val), [&](const Overflow_Map_Val& v){ return Overflow_Map_Val(v.cs, v.min_coord, f(v.val)); }, full);
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

    // Applies the overflow table's update with the k-mer `d`, holding no lock. Returns `false` iff it found the
    // overflow table full, having grown it instead; the update is to be retried then.
    const auto apply_overflow = [&](const Kmer_Desc& d) -> bool
    {
        bool full;
        const auto absent = ov_update(d.key, d.c, d.m, f, val, full);
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
            in_run = true, run_h = d.h, spilled = false;
            b = d.h & (cap_ - 1);
        }

        if(!spilled)
        {
            if(held == 0)
                lock(b), held = 1, cs_b = load_cs(md(b));

            uint32_t e;
            if(const auto j = probe(d.key, d.c, d.m, b, cs_b, e); j < B)
                return found(b, j), true;

            if(e < B)
            {
                place(d, b, e);
                cs_b = set_cs(cs_b, e, d.c);

                // A run's insertions fill its bucket's lines in order, and the slot prefetches predict the run's
                // first only: placing a key at a line's first slot prefetches the next line.
                if(e % line_slots == 0 && e + line_slots < B)
                    __builtin_prefetch(&slot(b, e + line_slots), 1, 3);

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
        if(const auto j = probe(d.key, d.c, d.m, b, e_b); j < B)
            return found(b, j), true;
        if(const auto j = probe(d.key, d.c, d.m, p, e_p); j < B)
            return found(p, j), true;
        if(const auto j = probe(d.key, d.c, d.m, q, e_q); j < B)
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
    stream<use_lanes, true, true, true>(seqs, lens, c,
        [&](const Kmer_Desc& d){ prefetch_meta<true>(d); },
        [&](const auto& make){ prefetch_slot<true>(make(), sc_slot); },
        [&](const Kmer_Desc& d, std::size_t, std::size_t)
        {
            if constexpr(mt_)
                if(KACHE_HASH_UNLIKELY(table_lock.writer_pending()))  // Yield to a resize.
                {
                    release();
                    table_lock.unlock_shared(token.id);
                    table_lock.lock_shared(token.id);
                    in_run = false;
                }

            while(!apply(d))
                ;
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
    bool sec = false;
    std::size_t b = 0, p = 0, q = 0;

    // Returns what `find` returns for the present k-mer at the `i`'th bucket's slot `j`.
    const auto present = [&](const std::size_t i, const std::size_t j) -> find_ret_t
    {
        (void)i, (void)j;
        if constexpr(is_map_)   return slot(i, j).second;
        else                    return true;
    };

    // Positive queries gain from prefetching their slots, and negative ones lose the work, so the slots are prefetched
    // while the last query hit; queries tend to come in streaks. A call's stretches go through a pipeline each, one
    // pipeline for all costs query performance. Lane runs don't continue across chunks for queries, which take no
    // locks: that costs them. Queries describe no minimizer coordinates: they don't use them.
    bool hit = true;
    Sec_Cache sc_slot;
    stream<use_lanes, false, false, false>(seqs, lens, c,
        [&](const Kmer_Desc& d){ prefetch_meta<false>(d); },
        Lookup_Pre_Slot{*this, sc_slot, hit},
        [&](const Kmer_Desc& d, const std::size_t seq_i, const std::size_t i)
        {
            // NB: Branch-free, as a branch on run starts mispredicts once per run.
            b = d.h & (cap_ - 1);
            sec &= (d.h == run_h);
            run_h = d.h;

            // Absent past a full primary bucket too, unless its spill filter marks the minimizer.
            uint32_t e;
            if(const auto j = probe(d.key, d.c, d.m, b, e); j < B)
                return hit = true, g(seq_i, i, present(b, j));
            if(e < B || !may_spill(md(b), d.h))
                return hit = false, g(seq_i, i, find_ret_t(null_val));

            if(!sec)
                std::tie(p, q) = secondaries(d.h), sec = true;

            uint32_t e_p, e_q;
            if(const auto j = probe(d.key, d.c, d.m, p, e_p); j < B)
                return hit = true, g(seq_i, i, present(p, j));
            if(const auto j = probe(d.key, d.c, d.m, q, e_q); j < B)
                return hit = true, g(seq_i, i, present(q, j));

            const auto r = (e_p == B && e_q == B) ? find_in_overflow(d.key) : find_ret_t(null_val);
            hit = static_cast<bool>(r);
            g(seq_i, i, r);
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
