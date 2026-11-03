//! Rust bindings to kache-hash's C interface: concurrent tables of k-mers, i.e. sets of k-mers, and maps from k-mers
//! to unsigned values of 8, 16, 32 or 64 bits. Non-ACGT characters split sequences, and lowercase bases count as
//! bases. A table keeps each k-mer in its canonical form, so that lookups find a k-mer in either orientation; oriented
//! inserts and masked counts tell a k-mer's occurrences as is from those as its reverse complement, e.g. for
//! forward-strand membership.
//!
//! A table is built through a [`Builder`], from any number of threads, each inserting through its own [`Inserter`];
//! then frozen into a [`Table`], which answers queries from any number of threads. A map's [`Updater`] updates the
//! values of the k-mers present, from any number of threads, while it borrows the table. The phases keep inserts,
//! updates and queries from overlapping, as the table requires.

#![deny(unsafe_op_in_unsafe_fn)]

use std::any::Any;
use std::ffi::CStr;
use std::fmt;
use std::marker::PhantomData;
use std::os::raw::{c_char, c_int, c_void};
use std::panic::{self, AssertUnwindSafe};
use std::ptr::NonNull;
use std::sync::Mutex;

/// Statistics of a table.
#[repr(C)]
#[derive(Debug, Clone, Copy, Default)]
pub struct Stats {
    /// Number of slots of the table.
    pub capacity: u64,
    /// Bytes of the table's buckets.
    pub bytes: u64,
    /// Number of keys.
    pub keys: u64,
    /// Number of keys in the overflow table.
    pub overflow_keys: u64,
    /// Percentages of the keys in their primary buckets, secondary buckets, and the overflow table.
    pub primary_pct: f64,
    pub secondary_pct: f64,
    pub overflow_pct: f64,
}

#[repr(C)]
#[derive(Clone, Copy)]
struct RawOptions {
    k: u32,
    value_bytes: u32,
    expected_keys: u64,
    load_factor: f64,
    resize_workers: u32,
    pow2_capacity: u32,
}

#[repr(C)]
struct RawTable {
    _private: [u8; 0],
}

type KeysFn = unsafe extern "C" fn(ctx: *mut c_void, keys: *const u64, values: *const c_void, n: usize) -> c_int;

extern "C" {
    fn kache_supports(k: u32, value_bytes: u32) -> c_int;
    fn kache_default_options(k: u32, value_bytes: u32) -> RawOptions;
    fn kache_new(options: *const RawOptions, out: *mut *mut RawTable) -> c_int;
    fn kache_free(t: *mut RawTable);
    fn kache_register(t: *mut RawTable, user: *mut u32) -> c_int;
    fn kache_insert(t: *mut RawTable, user: u32, seqs: *const *const c_char, lens: *const usize, n: usize, op: u32, x: u64, inserted: *mut u64) -> c_int;
    fn kache_insert_oriented(t: *mut RawTable, user: u32, seqs: *const *const c_char, lens: *const usize, n: usize, op: u32, x_fwd: u64, x_rev: u64, inserted: *mut u64) -> c_int;
    fn kache_update(t: *mut RawTable, seqs: *const *const c_char, lens: *const usize, n: usize, op: u32, x: u64, present: *mut u64) -> c_int;
    fn kache_count(t: *const RawTable, seqs: *const *const c_char, lens: *const usize, n: usize, counts: *mut u64) -> c_int;
    fn kache_count_masked(t: *const RawTable, seqs: *const *const c_char, lens: *const usize, n: usize, mask_fwd: u64, mask_rev: u64, counts: *mut u64) -> c_int;
    fn kache_sum(t: *const RawTable, seqs: *const *const c_char, lens: *const usize, n: usize, sums: *mut u64) -> c_int;
    fn kache_lookup(t: *const RawTable, seqs: *const *const c_char, lens: *const usize, n: usize, found: *mut u8, values: *mut c_void, fwd: *mut u8) -> c_int;
    fn kache_get(t: *const RawTable, kmer: *const c_char, found: *mut u8, value: *mut u64) -> c_int;
    fn kache_size(t: *const RawTable, keys: *mut u64) -> c_int;
    fn kache_reduce(t: *const RawTable, what: u32, mask: u64, workers: u32, out: *mut u64) -> c_int;
    fn kache_for_each(t: *const RawTable, shards: u32, shard: u32, f: KeysFn, ctx: *mut c_void) -> c_int;
    fn kache_get_stats(t: *const RawTable, out: *mut Stats) -> c_int;
    fn kache_strerror(code: c_int) -> *const c_char;
}

/// An error status of the C interface.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Error(pub i32);

impl Error {
    /// An invalid argument.
    pub const INVALID: Error = Error(1);
    /// No more users can register: a table takes 128 inserters at once.
    pub const USERS: Error = Error(2);
    /// Out of memory.
    pub const NO_MEMORY: Error = Error(3);
    /// Any other error of the table.
    pub const INTERNAL: Error = Error(4);
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // SAFETY: `kache_strerror` returns a static C string for any code.
        let msg = unsafe { CStr::from_ptr(kache_strerror(self.0)) };
        write!(f, "kache-hash: {} ({})", msg.to_string_lossy(), self.0)
    }
}

impl std::error::Error for Error {}

fn status(code: c_int) -> Result<(), Error> {
    if code == 0 {
        Ok(())
    } else {
        Err(Error(code))
    }
}

/// How an insert or an update combines a k-mer's value `v` with an operand `x`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u32)]
pub enum Op {
    /// `v`: a k-mer keeps the value it was inserted with.
    Keep = 0,
    /// `x`.
    Assign = 1,
    /// `v + x`, saturating at the values' maximum.
    Add = 2,
    /// `v | x`.
    Or = 3,
    /// `min(v, x)`.
    Min = 4,
    /// `max(v, x)`.
    Max = 5,
}

/// What a reduction over a map's keys computes, with a mask `m` of their values `v`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u32)]
pub enum Reduction {
    /// The number of keys with `v & m` nonzero.
    Keys = 0,
    /// The sum of the numbers of bits set in `v & m`.
    Popcount = 1,
    /// The sum of `v & m`, modulo 2^64.
    Sum = 2,
}

/// A lookup's finding at a k-mer position: the k-mer is absent.
pub const ABSENT: u8 = 0;
/// A lookup's finding at a k-mer position: the k-mer is present.
pub const PRESENT: u8 = 1;
/// A lookup's finding at a k-mer position: the position's k characters have a non-ACGT one.
pub const NOT_KMER: u8 = 2;

mod sealed {
    pub trait Sealed {}
}

/// The values of a table: `()` in a set, and `u8`, `u16`, `u32` or `u64` in a map.
pub trait Value: sealed::Sealed + Copy + Send + Sync + 'static {
    #[doc(hidden)]
    const BYTES: u32;
}

/// The values of a map.
pub trait MapValue: Value + Into<u64> {
    #[doc(hidden)]
    fn from_bits(x: u64) -> Self;
}

impl sealed::Sealed for () {}
impl Value for () {
    const BYTES: u32 = 0;
}

macro_rules! map_values {
    ($($t:ty),*) => {$(
        impl sealed::Sealed for $t {}
        impl Value for $t {
            const BYTES: u32 = std::mem::size_of::<$t>() as u32;
        }
        impl MapValue for $t {
            fn from_bits(x: u64) -> Self {
                x as $t
            }
        }
    )*};
}
map_values!(u8, u16, u32, u64);

/// Returns whether tables of `k`-mers with values of type `V` are supported: of 31-mers only, for now.
pub fn supports<V: Value>(k: usize) -> bool {
    // SAFETY: a pure function of its arguments.
    u32::try_from(k).is_ok_and(|k| unsafe { kache_supports(k, V::BYTES) } != 0)
}

/// Options of a table.
#[derive(Debug, Clone, Copy)]
pub struct Options {
    /// Length of the k-mers.
    pub k: usize,
    /// Number of k-mers to make room for; the table grows beyond, as needed.
    pub expected_keys: u64,
    /// Load factor, in `(0, 1)`, beyond which the table grows.
    pub load_factor: f64,
    /// Number of threads for the table's resizes and statistics, at least 1.
    pub resize_workers: usize,
    /// Whether to round the table's number of buckets up to a power of 2, for a lower load; else it is sized exactly
    /// for `expected_keys` at `load_factor`.
    pub pow2_capacity: bool,
}

impl Options {
    /// Returns the default options of a table of `k`-mers: room for 4M k-mers at load factor 0.8, its number of buckets
    /// rounded up to a power of 2, and 16 threads for resizes.
    pub fn new(k: usize) -> Self {
        // SAFETY: a pure function of its arguments.
        let o = unsafe { kache_default_options(u32::try_from(k).unwrap_or(0), 0) };
        Self {
            k,
            expected_keys: o.expected_keys,
            load_factor: o.load_factor,
            resize_workers: o.resize_workers as usize,
            pow2_capacity: o.pow2_capacity != 0,
        }
    }
}

/// Most sequences per call into the C interface: their pointers and lengths go on the stack.
const BATCH: usize = 64;

/// Calls `f` on the pointers and lengths of consecutive batches of at most `BATCH` of the sequences `seqs`, with the
/// index of each batch's first sequence.
fn batched<S: AsRef<[u8]>>(seqs: &[S], mut f: impl FnMut(&[*const c_char], &[usize], usize) -> Result<(), Error>) -> Result<(), Error> {
    let mut ptrs = [std::ptr::null::<c_char>(); BATCH];
    let mut lens = [0usize; BATCH];
    for (b, chunk) in seqs.chunks(BATCH).enumerate() {
        for (j, s) in chunk.iter().enumerate() {
            ptrs[j] = s.as_ref().as_ptr().cast();
            lens[j] = s.as_ref().len();
        }
        f(&ptrs[..chunk.len()], &lens[..chunk.len()], b * BATCH)?;
    }
    Ok(())
}

/// A table's handle, and its users that inserters have released, for new inserters.
struct Raw {
    t: NonNull<RawTable>,
    k: usize,
    free: Mutex<Vec<u32>>,
}

impl Raw {
    /// Returns the number of k-mer positions of the sequence `s`.
    fn positions(&self, s: &[u8]) -> usize {
        (s.len() + 1).saturating_sub(self.k)
    }
}

impl Drop for Raw {
    fn drop(&mut self) {
        // SAFETY: the table is freed once, after every borrow of it.
        unsafe { kache_free(self.t.as_ptr()) }
    }
}

/// A table being built: shared by the inserting threads, each through its own [`Inserter`].
pub struct Builder<V: Value> {
    raw: Raw,
    _v: PhantomData<V>,
}

// SAFETY: the C interface's inserts are thread-safe, each by its own user, and registrations are locked.
unsafe impl<V: Value> Send for Builder<V> {}
unsafe impl<V: Value> Sync for Builder<V> {}

impl<V: Value> Builder<V> {
    /// Creates a table with the options `options`.
    pub fn new(options: &Options) -> Result<Self, Error> {
        let o = RawOptions {
            k: u32::try_from(options.k).map_err(|_| Error::INVALID)?,
            value_bytes: V::BYTES,
            expected_keys: options.expected_keys,
            load_factor: options.load_factor,
            resize_workers: u32::try_from(options.resize_workers).map_err(|_| Error::INVALID)?,
            pow2_capacity: u32::from(options.pow2_capacity),
        };
        let mut t = std::ptr::null_mut();
        // SAFETY: `o` and `t` are valid for the call.
        status(unsafe { kache_new(&o, &mut t) })?;
        let t = NonNull::new(t).ok_or(Error::INTERNAL)?;
        Ok(Self { raw: Raw { t, k: options.k, free: Mutex::new(Vec::new()) }, _v: PhantomData })
    }

    /// Returns a new inserter: one per thread, and at most 128 at once.
    pub fn inserter(&self) -> Result<Inserter<'_, V>, Error> {
        let released = self.raw.free.lock().unwrap_or_else(|e| e.into_inner()).pop();
        let user = match released {
            Some(u) => u,
            None => {
                let mut u = 0;
                // SAFETY: the table is live, and registrations are locked.
                status(unsafe { kache_register(self.raw.t.as_ptr(), &mut u) })?;
                u
            }
        };
        Ok(Inserter { builder: self, user })
    }

    /// Freezes the table for queries, once every inserter is gone.
    pub fn freeze(self) -> Table<V> {
        Table { raw: self.raw, _v: PhantomData }
    }
}

/// An inserter of a table: one per thread.
pub struct Inserter<'a, V: Value> {
    builder: &'a Builder<V>,
    user: u32,
}

impl<V: Value> Inserter<'_, V> {
    fn insert_raw<S: AsRef<[u8]>>(&mut self, seqs: &[S], oriented: bool, op: u32, x_fwd: u64, x_rev: u64) -> Result<u64, Error> {
        let (t, user) = (self.builder.raw.t, self.user);
        let mut inserted = 0;
        batched(seqs, |p, l, _| {
            let mut n = 0;
            // SAFETY: the table is live while its builder is, the sequences are valid for their lengths, and this
            // inserter's user inserts from one thread at a time, as it is not shared.
            status(unsafe {
                if oriented {
                    kache_insert_oriented(t.as_ptr(), user, p.as_ptr(), l.as_ptr(), p.len(), op, x_fwd, x_rev, &mut n)
                } else {
                    kache_insert(t.as_ptr(), user, p.as_ptr(), l.as_ptr(), p.len(), op, x_fwd, &mut n)
                }
            })?;
            inserted += n;
            Ok(())
        })?;
        Ok(inserted)
    }
}

impl Inserter<'_, ()> {
    /// Inserts the k-mers of the sequences `seqs` into the set. Returns the number of k-mers newly inserted.
    pub fn insert<S: AsRef<[u8]>>(&mut self, seqs: &[S]) -> Result<u64, Error> {
        self.insert_raw(seqs, false, Op::Keep as u32, 0, 0)
    }
}

impl<V: MapValue> Inserter<'_, V> {
    /// Inserts the k-mers of the sequences `seqs` into the map: an absent k-mer gets the value `x`, and a present one's
    /// value `v` becomes `op(v, x)`. Returns the number of k-mers newly inserted.
    pub fn insert<S: AsRef<[u8]>>(&mut self, seqs: &[S], op: Op, x: V) -> Result<u64, Error> {
        self.insert_raw(seqs, false, op as u32, x.into(), 0)
    }

    /// Inserts the k-mers of the sequences `seqs` into the map as [`Inserter::insert`] does, with the operand `x_fwd`
    /// for each occurrence of a k-mer as is, and `x_rev` for each as its reverse complement. The operation `op` is one
    /// of [`Op::Assign`], [`Op::Add`], [`Op::Or`] and [`Op::Max`].
    pub fn insert_oriented<S: AsRef<[u8]>>(&mut self, seqs: &[S], op: Op, x_fwd: V, x_rev: V) -> Result<u64, Error> {
        self.insert_raw(seqs, true, op as u32, x_fwd.into(), x_rev.into())
    }
}

impl<V: Value> Drop for Inserter<'_, V> {
    fn drop(&mut self) {
        // The user's inserts are over; another inserter may take it, on any thread.
        self.builder.raw.free.lock().unwrap_or_else(|e| e.into_inner()).push(self.user);
    }
}

/// A frozen table: queried from any number of threads.
pub struct Table<V: Value> {
    raw: Raw,
    _v: PhantomData<V>,
}

// SAFETY: the C interface's queries are thread-safe, and nothing inserts into a frozen table, nor updates it but
// through an `Updater`, which borrows it mutably.
unsafe impl<V: Value> Send for Table<V> {}
unsafe impl<V: Value> Sync for Table<V> {}

/// The context of an iteration's callbacks: the receiver, and the panic it raised, if any.
struct Chunks<'f, V, F> {
    f: &'f mut F,
    panic: Option<Box<dyn Any + Send>>,
    _v: PhantomData<V>,
}

/// Passes a chunk of keys and values to the receiver of the iteration `ctx`; stops the iteration if it asks to, or if
/// it panics, keeping the panic for the caller.
unsafe extern "C" fn receive<V: Value, F: FnMut(&[u64], &[V]) -> bool>(ctx: *mut c_void, keys: *const u64, values: *const c_void, n: usize) -> c_int {
    // SAFETY: `ctx` is the iteration's `Chunks`, alive and unshared during the call, and the C interface passes `n`
    // keys, and `n` values of `V`'s width in a map; the values of a set are zero-sized.
    let c = unsafe { &mut *ctx.cast::<Chunks<'_, V, F>>() };
    let keys = unsafe { std::slice::from_raw_parts(keys, n) };
    let values: &[V] = if V::BYTES == 0 {
        unsafe { std::slice::from_raw_parts(NonNull::<V>::dangling().as_ptr(), n) }
    } else {
        unsafe { std::slice::from_raw_parts(values.cast::<V>(), n) }
    };
    match panic::catch_unwind(AssertUnwindSafe(|| (c.f)(keys, values))) {
        Ok(go_on) => c_int::from(!go_on),
        Err(p) => {
            c.panic = Some(p);
            1
        }
    }
}

impl<V: Value> Table<V> {
    /// Returns the length of the table's k-mers.
    pub fn k(&self) -> usize {
        self.raw.k
    }

    /// Returns the number of k-mer positions of the sequences `seqs`: the size of their lookups' outputs.
    pub fn positions<S: AsRef<[u8]>>(&self, seqs: &[S]) -> usize {
        seqs.iter().map(|s| self.raw.positions(s.as_ref())).sum()
    }

    /// Puts the number of k-mers of each sequence of `seqs` present in the table at the matching index of `counts`:
    /// k-mer positions, with repeats counted.
    pub fn count<S: AsRef<[u8]>>(&self, seqs: &[S], counts: &mut [u64]) -> Result<(), Error> {
        assert!(counts.len() >= seqs.len(), "a count for each sequence");
        let t = self.raw.t;
        batched(seqs, |p, l, i| {
            // SAFETY: the table is live, the sequences are valid for their lengths, and `counts[i..]` has room for the
            // batch.
            status(unsafe { kache_count(t.as_ptr(), p.as_ptr(), l.as_ptr(), p.len(), counts[i..].as_mut_ptr()) })
        })
    }

    /// Looks up each k-mer position of the sequences `seqs`, in order: puts [`PRESENT`], [`ABSENT`] or [`NOT_KMER`] in
    /// `found`, and whether the k-mer occurs as is in `fwd`, if given. Each has room for [`Table::positions`].
    pub fn lookup<S: AsRef<[u8]>>(&self, seqs: &[S], found: &mut [u8], fwd: Option<&mut [u8]>) -> Result<(), Error> {
        self.lookup_raw(seqs, found, std::ptr::null_mut(), 0, fwd)
    }

    fn lookup_raw<S: AsRef<[u8]>>(&self, seqs: &[S], found: &mut [u8], values: *mut V, values_len: usize, mut fwd: Option<&mut [u8]>) -> Result<(), Error> {
        let p = self.positions(seqs);
        assert!(found.len() >= p && (values.is_null() || values_len >= p), "an output for each k-mer position");
        assert!(fwd.as_ref().is_none_or(|f| f.len() >= p), "an orientation for each k-mer position");
        let t = self.raw.t;
        let mut o = 0;  // Index of the batch's first position in the outputs.
        batched(seqs, |ptrs, lens, _| {
            let v = if values.is_null() { std::ptr::null_mut() } else { values.wrapping_add(o).cast::<c_void>() };
            let f = fwd.as_mut().map_or(std::ptr::null_mut(), |f| f[o..].as_mut_ptr());
            // SAFETY: the table is live, the sequences are valid for their lengths, and the outputs have room for the
            // batch's positions from `o`.
            status(unsafe { kache_lookup(t.as_ptr(), ptrs.as_ptr(), lens.as_ptr(), ptrs.len(), found[o..].as_mut_ptr(), v, f) })?;
            o += lens.iter().map(|&n| (n + 1).saturating_sub(self.raw.k)).sum::<usize>();
            Ok(())
        })
    }

    /// Returns whether the k-mer of the first `k` characters of `kmer` is present, in either orientation; one with a
    /// non-ACGT character is absent.
    pub fn contains(&self, kmer: &[u8]) -> Result<bool, Error> {
        self.get_raw(kmer).map(|(found, _)| found)
    }

    fn get_raw(&self, kmer: &[u8]) -> Result<(bool, u64), Error> {
        assert!(kmer.len() >= self.raw.k, "a k-mer of k characters");
        let (mut found, mut v) = (0, 0);
        // SAFETY: the table is live, and `kmer` has at least `k` characters.
        status(unsafe { kache_get(self.raw.t.as_ptr(), kmer.as_ptr().cast(), &mut found, &mut v) })?;
        Ok((found != 0, v))
    }

    /// Returns the number of keys of the table, i.e. of distinct canonical k-mers.
    pub fn len(&self) -> u64 {
        let mut n = 0;
        // SAFETY: the table is live, and `n` is a valid out-pointer.
        status(unsafe { kache_size(self.raw.t.as_ptr(), &mut n) }).expect("the size of a live table");
        n
    }

    /// Returns whether the table has no keys.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Calls `f(keys, values)` with the keys of the table in the share `shard` of `shards` shares, in chunks and in an
    /// unspecified order, while it returns `true`: each key in its 2-bit encoding, with A, C, G and T as 0 to 3 and the
    /// first base the most significant, and in a map, with its value. The shares partition the table, so that
    /// iterations over all of them, concurrent or not, visit each key once. A panic of `f` ends the iteration, and
    /// resumes in the caller.
    pub fn for_each_chunk<F: FnMut(&[u64], &[V]) -> bool>(&self, shards: usize, shard: usize, mut f: F) -> Result<(), Error> {
        let shards = u32::try_from(shards).map_err(|_| Error::INVALID)?;
        let shard = u32::try_from(shard).map_err(|_| Error::INVALID)?;
        let mut c = Chunks::<V, F> { f: &mut f, panic: None, _v: PhantomData };
        // SAFETY: the table is live, and `c` outlives the iteration, which calls `receive` with it alone.
        let r = status(unsafe { kache_for_each(self.raw.t.as_ptr(), shards, shard, receive::<V, F>, (&mut c as *mut Chunks<'_, V, F>).cast()) });
        if let Some(p) = c.panic {
            panic::resume_unwind(p);
        }
        r
    }

    /// Returns the statistics of the table.
    pub fn stats(&self) -> Result<Stats, Error> {
        let mut s = Stats::default();
        // SAFETY: the table is live, and `s` is a valid out-pointer.
        status(unsafe { kache_get_stats(self.raw.t.as_ptr(), &mut s) })?;
        Ok(s)
    }

    /// Thaws the table for more inserts.
    pub fn thaw(self) -> Builder<V> {
        Builder { raw: self.raw, _v: PhantomData }
    }
}

impl<V: MapValue> Table<V> {
    /// Puts the number of k-mers of each sequence of `seqs` present in the map with values matching their
    /// orientation's mask at the matching index of `counts`: of the k-mers with values `v`, the occurrences as is with
    /// `v & mask_fwd` nonzero, and those as reverse complements with `v & mask_rev` nonzero.
    pub fn count_masked<S: AsRef<[u8]>>(&self, seqs: &[S], mask_fwd: V, mask_rev: V, counts: &mut [u64]) -> Result<(), Error> {
        assert!(counts.len() >= seqs.len(), "a count for each sequence");
        let (t, m_fwd, m_rev) = (self.raw.t, mask_fwd.into(), mask_rev.into());
        batched(seqs, |p, l, i| {
            // SAFETY: as in `count`.
            status(unsafe { kache_count_masked(t.as_ptr(), p.as_ptr(), l.as_ptr(), p.len(), m_fwd, m_rev, counts[i..].as_mut_ptr()) })
        })
    }

    /// Puts the sum of the values of the k-mers of each sequence of `seqs` present in the map at the matching index of
    /// `sums`, with repeats counted, modulo 2^64.
    pub fn sum<S: AsRef<[u8]>>(&self, seqs: &[S], sums: &mut [u64]) -> Result<(), Error> {
        assert!(sums.len() >= seqs.len(), "a sum for each sequence");
        let t = self.raw.t;
        batched(seqs, |p, l, i| {
            // SAFETY: as in `count`.
            status(unsafe { kache_sum(t.as_ptr(), p.as_ptr(), l.as_ptr(), p.len(), sums[i..].as_mut_ptr()) })
        })
    }

    /// Looks up each k-mer position of the sequences `seqs` as [`Table::lookup`] does, and puts the k-mers' values, or
    /// 0, in `values`, which has room for [`Table::positions`] too.
    pub fn lookup_values<S: AsRef<[u8]>>(&self, seqs: &[S], found: &mut [u8], values: &mut [V], fwd: Option<&mut [u8]>) -> Result<(), Error> {
        self.lookup_raw(seqs, found, values.as_mut_ptr(), values.len(), fwd)
    }

    /// Returns the value of the k-mer of the first `k` characters of `kmer`, in either orientation, if present; one
    /// with a non-ACGT character is absent.
    pub fn get(&self, kmer: &[u8]) -> Result<Option<V>, Error> {
        let (found, v) = self.get_raw(kmer)?;
        Ok(found.then(|| V::from_bits(v)))
    }

    /// Returns the reduction `what` of the map's keys, with the mask `mask` of their values, computed by `workers`
    /// threads, at least 1.
    pub fn reduce(&self, what: Reduction, mask: V, workers: usize) -> Result<u64, Error> {
        let workers = u32::try_from(workers).map_err(|_| Error::INVALID)?;
        let mut r = 0;
        // SAFETY: the table is live, and `r` is a valid out-pointer.
        status(unsafe { kache_reduce(self.raw.t.as_ptr(), what as u32, mask.into(), workers, &mut r) })?;
        Ok(r)
    }

    /// Returns an updater of the map, which borrows it until dropped.
    pub fn updater(&mut self) -> Updater<'_, V> {
        Updater { table: self }
    }
}

/// An updater of a map: it updates the values of the k-mers present, from any number of threads at once, and keeps
/// the map from being queried meanwhile.
pub struct Updater<'a, V: MapValue> {
    table: &'a mut Table<V>,
}

impl<V: MapValue> Updater<'_, V> {
    /// Updates the value `v` of each k-mer of the sequences `seqs` present in the map to `op(v, x)`, and skips the
    /// absent ones: it inserts none. Returns the number of k-mers present, with repeats counted.
    pub fn update<S: AsRef<[u8]>>(&self, seqs: &[S], op: Op, x: V) -> Result<u64, Error> {
        let (t, x) = (self.table.raw.t, x.into());
        let mut present = 0;
        batched(seqs, |p, l, _| {
            let mut n = 0;
            // SAFETY: the map is live and borrowed by this updater, so that only updates run on it, which may run
            // concurrently; the sequences are valid for their lengths.
            status(unsafe { kache_update(t.as_ptr(), p.as_ptr(), l.as_ptr(), p.len(), op as u32, x, &mut n) })?;
            present += n;
            Ok(())
        })?;
        Ok(present)
    }
}
