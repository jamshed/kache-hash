
#ifndef KACHE_HASH_CAPI_H
#define KACHE_HASH_CAPI_H



// A C interface to kache-hash's concurrent tables of k-mers, for callers in
// other languages, such as Rust: a set of k-mers, or a map from k-mers to
// unsigned values of 8, 16, 32 or 64 bits. The operations take batches of
// sequences: non-ACGT characters split them, and lowercase bases count as
// bases. A table keeps each k-mer in its canonical form, the lesser of it and
// its reverse complement, so that lookups find a k-mer in either orientation.
// The oriented operations tell the orientations apart: a k-mer occurs as is,
// i.e. as its canonical form, or as that form's reverse complement.
//
// Inserts may run concurrently, each user on its own thread. Updates may run
// concurrently with one another, but with no other operations. Queries, i.e.
// counts, sums, lookups, reductions, iterations, sizes and statistics, may run
// concurrently with one another, but not with inserts or updates. The
// functions that can fail return a status code, `KACHE_OK` on success, and
// catch every exception of the table.

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#define KACHE_NOEXCEPT noexcept
#else
#define KACHE_NOEXCEPT
#endif


enum
{
    KACHE_OK = 0,
    KACHE_EINVAL = 1,       // An invalid argument.
    KACHE_EUSERS = 2,       // No more users can register.
    KACHE_ENOMEM = 3,       // Out of memory.
    KACHE_EINTERNAL = 4     // Any other error of the table.
};

#define KACHE_MAX_USERS 128 // Most users a table registers.

// How an insert or an update combines a k-mer's value `v` with an operand `x`.
enum
{
    KACHE_KEEP = 0,     // `v`: a k-mer keeps the value it was inserted with.
    KACHE_ASSIGN = 1,   // `x`.
    KACHE_ADD = 2,      // `v + x`, saturating at the values' maximum.
    KACHE_OR = 3,       // `v | x`.
    KACHE_MIN = 4,      // `min(v, x)`.
    KACHE_MAX = 5       // `max(v, x)`.
};

// What a reduction over a table's keys computes, with a mask `m` of their
// values `v`.
enum
{
    KACHE_REDUCE_KEYS = 0,      // The number of keys with `v & m` nonzero; of all keys in a set.
    KACHE_REDUCE_POPCOUNT = 1,  // The sum of the numbers of bits set in `v & m`.
    KACHE_REDUCE_SUM = 2        // The sum of `v & m`, modulo 2^64.
};

// What a lookup finds at a k-mer position.
enum
{
    KACHE_ABSENT = 0,
    KACHE_PRESENT = 1,
    KACHE_NOT_KMER = 2  // The position's k characters have a non-ACGT one.
};


// A table: a set of k-mers, or a map from k-mers to values.
typedef struct kache_table kache_table;

// Options of a table.
typedef struct kache_options
{
    uint32_t k;                 // Length of the k-mers.
    uint32_t value_bytes;       // Width of the values in bytes: 1, 2, 4 or 8 in a map, and 0 in a set.
    uint64_t expected_keys;     // Number of k-mers to make room for; the table grows beyond, as needed.
    double load_factor;         // Load factor, in `(0, 1)`, beyond which the table grows.
    uint32_t resize_workers;    // Number of threads for the table's resizes and statistics, at least 1.
    uint32_t pow2_capacity;     // 1 to round the table's number of buckets up to a power of 2, for a lower load; 0 to
                                // size it exactly for `expected_keys` at `load_factor`.
} kache_options;

// Statistics of a table.
typedef struct kache_stats
{
    uint64_t capacity;      // Number of slots of the table.
    uint64_t bytes;         // Bytes of the table's buckets.
    uint64_t keys;          // Number of keys.
    uint64_t overflow_keys; // Number of keys in the overflow table.
    double primary_pct;     // Percentages of the keys in their primary buckets, secondary buckets, and the overflow table.
    double secondary_pct;
    double overflow_pct;
} kache_stats;

// Receives a chunk of a table's keys in an iteration: `keys[0..n)`, each in
// its 2-bit encoding, with A, C, G and T as 0 to 3 and the first base the most
// significant, and in a map, their values `values[0..n)` of the table's width;
// NULL in a set. Returns 0 to go on, or nonzero to receive no more chunks. It
// must not modify the table.
typedef int (*kache_keys_fn)(void* ctx, const uint64_t* keys, const void* values, size_t n);


// Returns whether tables of `k`-mers with values of `value_bytes` bytes are
// supported: of 31-mers only, for now, with each width.
int kache_supports(uint32_t k, uint32_t value_bytes) KACHE_NOEXCEPT;

// Returns the default options of a table of `k`-mers with values of
// `value_bytes` bytes: room for 4M k-mers at load factor 0.8, sized exactly,
// and 16 threads for resizes.
kache_options kache_default_options(uint32_t k, uint32_t value_bytes) KACHE_NOEXCEPT;

// Creates a table with the options `*options` at `*out`. A table takes at most
// 2^32 - 1 buckets of 32 k-mers, or 2^31 with `pow2_capacity`.
int kache_new(const kache_options* options, kache_table** out) KACHE_NOEXCEPT;

// Frees the table `t`; nothing if `t` is NULL.
void kache_free(kache_table* t) KACHE_NOEXCEPT;

// Registers a new user of the table `t` for inserts at `*user`: at most
// `KACHE_MAX_USERS` per table. A user inserts from one thread at a time, and
// may pass to another thread between its inserts.
int kache_register(kache_table* t, uint32_t* user) KACHE_NOEXCEPT;

// Inserts the k-mers of the sequences `seqs[0..n)` of lengths `lens[0..n)` into
// the table `t`, by the registered user `user`. In a map, an absent k-mer gets
// the value `x`, and a present one's value `v` becomes `op(v, x)`, by the
// operation `op`, one of `KACHE_KEEP` to `KACHE_MAX`; a set ignores both. Puts
// the number of k-mers inserted at `*inserted`, unless it is NULL.
int kache_insert(kache_table* t, uint32_t user, const char* const* seqs, const size_t* lens, size_t n, uint32_t op, uint64_t x, uint64_t* inserted) KACHE_NOEXCEPT;

// Inserts the k-mers of the sequences into the map `t` as `kache_insert` does,
// with the operand `x_fwd` for each occurrence of a k-mer as is, and `x_rev`
// for each as its reverse complement. The operation `op` is one of
// `KACHE_ASSIGN`, `KACHE_ADD`, `KACHE_OR` and `KACHE_MAX`.
int kache_insert_oriented(kache_table* t, uint32_t user, const char* const* seqs, const size_t* lens, size_t n, uint32_t op, uint64_t x_fwd, uint64_t x_rev, uint64_t* inserted) KACHE_NOEXCEPT;

// Updates the value `v` of each k-mer of the sequences `seqs[0..n)` of lengths
// `lens[0..n)` present in the map `t` to `op(v, x)`, and skips the absent ones:
// it inserts none. Puts the number of k-mers present, with repeats counted, at
// `*present`, unless it is NULL. It takes no user, nor locks.
int kache_update(kache_table* t, const char* const* seqs, const size_t* lens, size_t n, uint32_t op, uint64_t x, uint64_t* present) KACHE_NOEXCEPT;

// Puts the number of k-mers of the sequence `seqs[i]`, of length `lens[i]`,
// present in the table `t` at `counts[i]`, for each `i` in `[0, n)`: k-mer
// positions, with repeats counted.
int kache_count(const kache_table* t, const char* const* seqs, const size_t* lens, size_t n, uint64_t* counts) KACHE_NOEXCEPT;

// Puts the number of k-mers of each sequence present in the map `t` whose
// values match their orientation's mask at `counts[i]`, as `kache_count`
// does: of the k-mers with values `v`, the occurrences as is with
// `v & mask_fwd` nonzero, and those as reverse complements with
// `v & mask_rev` nonzero.
int kache_count_masked(const kache_table* t, const char* const* seqs, const size_t* lens, size_t n, uint64_t mask_fwd, uint64_t mask_rev, uint64_t* counts) KACHE_NOEXCEPT;

// Puts the sum of the values of the k-mers of the sequence `seqs[i]`, of
// length `lens[i]`, present in the map `t` at `sums[i]`, for each `i` in
// `[0, n)`, with repeats counted, modulo 2^64.
int kache_sum(const kache_table* t, const char* const* seqs, const size_t* lens, size_t n, uint64_t* sums) KACHE_NOEXCEPT;

// Looks up each k-mer position of the sequences `seqs[0..n)` of lengths
// `lens[0..n)` in the table `t`. The position `j` of `seqs[i]` is at index
// `o_i + j` of the outputs, where `o_i` is the number of positions of the
// sequences before it: the sum of `lens[h] - k + 1` over the sequences `h < i`
// of at least `k` characters. Puts `KACHE_PRESENT`, `KACHE_ABSENT` or
// `KACHE_NOT_KMER` at `found[o_i + j]`; in a map, the k-mer's value, or 0, at
// `values[o_i + j]`, of the table's width, unless `values` is NULL; and
// whether the k-mer occurs as is at `fwd[o_i + j]`, unless `fwd` is NULL.
int kache_lookup(const kache_table* t, const char* const* seqs, const size_t* lens, size_t n, uint8_t* found, void* values, uint8_t* fwd) KACHE_NOEXCEPT;

// Looks up the k-mer of the `k` characters at `kmer` in the table `t`, in
// either orientation: puts whether it is present at `*found`, and in a map,
// its value, or 0, at `*value`, unless it is NULL. A k-mer with a non-ACGT
// character is absent.
int kache_get(const kache_table* t, const char* kmer, uint8_t* found, uint64_t* value) KACHE_NOEXCEPT;

// Puts the number of keys of the table `t`, i.e. of distinct canonical k-mers,
// at `*keys`.
int kache_size(const kache_table* t, uint64_t* keys) KACHE_NOEXCEPT;

// Puts the reduction `what` of the keys of the table `t`, one of
// `KACHE_REDUCE_KEYS` to `KACHE_REDUCE_SUM`, with the mask `mask` of their
// values, at `*out`, computed by `workers` threads, at least 1. A set supports
// `KACHE_REDUCE_KEYS` only, and ignores the mask.
int kache_reduce(const kache_table* t, uint32_t what, uint64_t mask, uint32_t workers, uint64_t* out) KACHE_NOEXCEPT;

// Calls `f(ctx, keys, values, n)` with the keys of the table `t` in the share
// `shard` of `shards` shares, in chunks, in an unspecified order, until it
// returns nonzero. The shares partition the table, so that iterations over
// all of them, concurrent or not, visit each key once.
int kache_for_each(const kache_table* t, uint32_t shards, uint32_t shard, kache_keys_fn f, void* ctx) KACHE_NOEXCEPT;

// Puts the statistics of the table `t` at `*out`.
int kache_get_stats(const kache_table* t, kache_stats* out) KACHE_NOEXCEPT;

// Returns a description of the status code `code`.
const char* kache_strerror(int code) KACHE_NOEXCEPT;


#ifdef __cplusplus
}
#endif



#endif
