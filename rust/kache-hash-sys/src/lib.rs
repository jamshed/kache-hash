//! Raw bindings to kache-hash's C interface, `kache_hash/kache_capi.h`, which documents each function, and the
//! interface library's build: concurrent tables of k-mers, i.e. sets of k-mers, and maps from k-mers to unsigned values
//! of 8, 16, 32 or 64 bits. The crate `kache-hash` wraps them in a safe interface.
//!
//! The build compiles the library with the `cc` crate, for an x86-64 CPU with AVX2 at least: for the ISA that
//! `KACHE_HASH_ARCH` names, if set, `native`, `x86-64-v3`, `x86-64-v4`, or `none` for the compiler's own flags;
//! otherwise for the native CPU if rustc targets it (`-C target-cpu=native`), for x86-64-v4 if rustc's target features
//! include AVX-512's, and for x86-64-v3 otherwise.

#![no_std]
#![allow(non_camel_case_types)]

use core::ffi::{c_char, c_int, c_void};
use core::marker::{PhantomData, PhantomPinned};

/// A status: success.
pub const KACHE_OK: c_int = 0;
/// A status: an invalid argument.
pub const KACHE_EINVAL: c_int = 1;
/// A status: no more users can register.
pub const KACHE_EUSERS: c_int = 2;
/// A status: out of memory.
pub const KACHE_ENOMEM: c_int = 3;
/// A status: any other error of the table.
pub const KACHE_EINTERNAL: c_int = 4;

/// Most users a table registers.
pub const KACHE_MAX_USERS: u32 = 128;

/// An operation of inserts and updates on a k-mer's value `v` and an operand `x`: `v`, which keeps the value that the
/// k-mer was inserted with.
pub const KACHE_KEEP: u32 = 0;
/// An operation: `x`.
pub const KACHE_ASSIGN: u32 = 1;
/// An operation: `v + x`, saturating at the values' maximum.
pub const KACHE_ADD: u32 = 2;
/// An operation: `v | x`.
pub const KACHE_OR: u32 = 3;
/// An operation: `min(v, x)`.
pub const KACHE_MIN: u32 = 4;
/// An operation: `max(v, x)`.
pub const KACHE_MAX: u32 = 5;

/// A reduction over a table's keys, with a mask `m` of their values `v`: the number of keys with `v & m` nonzero; of
/// all keys in a set.
pub const KACHE_REDUCE_KEYS: u32 = 0;
/// A reduction: the sum of the numbers of bits set in `v & m`.
pub const KACHE_REDUCE_POPCOUNT: u32 = 1;
/// A reduction: the sum of `v & m`, modulo 2^64.
pub const KACHE_REDUCE_SUM: u32 = 2;

/// A lookup's finding at a k-mer position: the k-mer is absent.
pub const KACHE_ABSENT: u8 = 0;
/// A lookup's finding: the k-mer is present.
pub const KACHE_PRESENT: u8 = 1;
/// A lookup's finding: the position's k characters have a non-ACGT one.
pub const KACHE_NOT_KMER: u8 = 2;

/// A table: a set of k-mers, or a map from k-mers to values.
#[repr(C)]
pub struct kache_table {
    _data: [u8; 0],
    _marker: PhantomData<(*mut u8, PhantomPinned)>,
}

/// Options of a table.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct kache_options {
    /// Length of the k-mers.
    pub k: u32,
    /// Width of the values in bytes: 1, 2, 4 or 8 in a map, and 0 in a set.
    pub value_bytes: u32,
    /// Number of k-mers to make room for; the table grows beyond, as needed.
    pub expected_keys: u64,
    /// Load factor, in `(0, 1)`, beyond which the table grows.
    pub load_factor: f64,
    /// Number of threads for the table's resizes and statistics, at least 1.
    pub resize_workers: u32,
    /// 0, the default, to size the table exactly for `expected_keys` at `load_factor`; 1 to round its number of
    /// buckets up to a power of 2, for a lower load.
    pub pow2_capacity: u32,
}

/// Statistics of a table.
#[repr(C)]
#[derive(Debug, Clone, Copy, Default)]
pub struct kache_stats {
    /// Number of slots of the table.
    pub capacity: u64,
    /// Bytes of the table's buckets.
    pub bytes: u64,
    /// Number of keys.
    pub keys: u64,
    /// Number of keys in the overflow table.
    pub overflow_keys: u64,
    /// Percentage of the keys in their primary buckets.
    pub primary_pct: f64,
    /// Percentage of the keys in their secondary buckets.
    pub secondary_pct: f64,
    /// Percentage of the keys in the overflow table.
    pub overflow_pct: f64,
}

/// The receiver of a chunk of a table's keys in an iteration, `kache_for_each`'s: `keys[0..n)`, in their 2-bit
/// encoding, and in a map their values, `values[0..n)`, of the table's width. It returns 0 to go on.
pub type kache_keys_fn =
    Option<unsafe extern "C" fn(ctx: *mut c_void, keys: *const u64, values: *const c_void, n: usize) -> c_int>;

extern "C" {
    pub fn kache_version() -> *const c_char;
    pub fn kache_supports(k: u32, value_bytes: u32) -> c_int;
    pub fn kache_default_options(k: u32, value_bytes: u32) -> kache_options;
    pub fn kache_new(options: *const kache_options, out: *mut *mut kache_table) -> c_int;
    pub fn kache_free(t: *mut kache_table);
    pub fn kache_register(t: *mut kache_table, user: *mut u32) -> c_int;
    pub fn kache_insert(t: *mut kache_table, user: u32, seqs: *const *const c_char, lens: *const usize, n: usize, op: u32, x: u64, inserted: *mut u64) -> c_int;
    pub fn kache_insert_oriented(t: *mut kache_table, user: u32, seqs: *const *const c_char, lens: *const usize, n: usize, op: u32, x_fwd: u64, x_rev: u64, inserted: *mut u64) -> c_int;
    pub fn kache_update(t: *mut kache_table, seqs: *const *const c_char, lens: *const usize, n: usize, op: u32, x: u64, present: *mut u64) -> c_int;
    pub fn kache_count(t: *const kache_table, seqs: *const *const c_char, lens: *const usize, n: usize, counts: *mut u64) -> c_int;
    pub fn kache_count_masked(t: *const kache_table, seqs: *const *const c_char, lens: *const usize, n: usize, mask_fwd: u64, mask_rev: u64, counts: *mut u64) -> c_int;
    pub fn kache_sum(t: *const kache_table, seqs: *const *const c_char, lens: *const usize, n: usize, sums: *mut u64) -> c_int;
    pub fn kache_lookup(t: *const kache_table, seqs: *const *const c_char, lens: *const usize, n: usize, found: *mut u8, values: *mut c_void, fwd: *mut u8) -> c_int;
    pub fn kache_get(t: *const kache_table, kmer: *const c_char, found: *mut u8, value: *mut u64) -> c_int;
    pub fn kache_size(t: *const kache_table, keys: *mut u64) -> c_int;
    pub fn kache_reduce(t: *const kache_table, what: u32, mask: u64, workers: u32, out: *mut u64) -> c_int;
    pub fn kache_for_each(t: *const kache_table, shards: u32, shard: u32, f: kache_keys_fn, ctx: *mut c_void) -> c_int;
    pub fn kache_get_stats(t: *const kache_table, out: *mut kache_stats) -> c_int;
    pub fn kache_strerror(code: c_int) -> *const c_char;
}
