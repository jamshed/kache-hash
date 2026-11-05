//! C interface to a CBL k-mer set, for the benchmark driver `bench-cbl.cpp`.

#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use cbl::CBL;
use std::slice;

const K: usize = 31; // k-mer length; must match the driver's.
type T = u128; // Integer type of the k-mers; must fit 2k + lg(2k) bits.
const PREFIX_BITS: usize = 28; // Bits of the k-mers' prefixes, as CBL's paper sets for large inputs.

type Set = CBL<K, T, PREFIX_BITS>;

/// Returns a new empty set of canonical k-mers.
#[no_mangle]
pub extern "C" fn cbl_new() -> *mut Set {
    Box::into_raw(Box::new(Set::new_canonical()))
}

/// Frees the set `s`.
#[no_mangle]
pub unsafe extern "C" fn cbl_free(s: *mut Set) {
    drop(Box::from_raw(s));
}

/// Adds the k-mers of the sequence `seq` of length `len` ≥ k to the set `s`.
#[no_mangle]
pub unsafe extern "C" fn cbl_insert_seq(s: *mut Set, seq: *const u8, len: usize) {
    (*s).insert_seq(slice::from_raw_parts(seq, len));
}

/// Returns the number of the k-mers of the sequence `seq` of length `len` ≥ k
/// present in the set `s`.
#[no_mangle]
pub unsafe extern "C" fn cbl_count_seq(s: *mut Set, seq: *const u8, len: usize) -> usize {
    (*s).contains_seq(slice::from_raw_parts(seq, len))
        .iter()
        .filter(|&&b| b)
        .count()
}

/// Returns the number of k-mers in the set `s`.
#[no_mangle]
pub unsafe extern "C" fn cbl_count(s: *const Set) -> usize {
    (*s).count()
}

/// Returns the proportion of the prefixes used in the set `s`.
#[no_mangle]
pub unsafe extern "C" fn cbl_prefix_load(s: *const Set) -> f64 {
    (*s).prefix_load()
}

/// Returns the k-mer length.
#[no_mangle]
pub extern "C" fn cbl_k() -> usize {
    K
}

/// Returns the bits of the k-mers' prefixes.
#[no_mangle]
pub extern "C" fn cbl_prefix_bits() -> usize {
    PREFIX_BITS
}
