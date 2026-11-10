// Checks that the bindings match the library they link: its version is the crate's, its default options and
// statistics come through as the header lays them out, and a set takes and counts a sequence's k-mers.

use kache_hash_sys::*;
use std::ffi::CStr;
use std::ptr;

#[test]
fn version() {
    // SAFETY: `kache_version` returns a static C string.
    let v = unsafe { CStr::from_ptr(kache_version()) };
    assert_eq!(v.to_str().unwrap(), env!("CARGO_PKG_VERSION"));
}

#[test]
fn set() {
    let seq = b"ACGTACGGTCAGTTAGCATCGATCGGATCCATGCAGTCAT";
    let (seqs, lens) = ([seq.as_ptr().cast()], [seq.len()]);
    // SAFETY: the table is live from `kache_new` to `kache_free`, and the pointers are valid for each call.
    unsafe {
        assert_ne!(kache_supports(31, 0), 0);
        let mut o = kache_default_options(31, 0);
        assert_eq!((o.k, o.value_bytes, o.expected_keys, o.load_factor, o.resize_workers, o.pow2_capacity), (31, 0, 1 << 22, 0.8, 16, 0));

        o.expected_keys = 1000;
        let mut t = ptr::null_mut();
        assert_eq!(kache_new(&o, &mut t), KACHE_OK);
        let mut user = 0;
        assert_eq!(kache_register(t, &mut user), KACHE_OK);
        let mut inserted = 0;
        assert_eq!(kache_insert(t, user, seqs.as_ptr(), lens.as_ptr(), 1, KACHE_KEEP, 0, &mut inserted), KACHE_OK);
        assert_eq!(inserted, 10);

        let mut counts = [0];
        assert_eq!(kache_count(t, seqs.as_ptr(), lens.as_ptr(), 1, counts.as_mut_ptr()), KACHE_OK);
        assert_eq!(counts, [10]);
        let mut s = kache_stats::default();
        assert_eq!(kache_get_stats(t, &mut s), KACHE_OK);
        assert_eq!(s.keys, 10);
        assert!(s.capacity >= 1000 && s.bytes > 0);
        assert!((s.primary_pct + s.secondary_pct + s.overflow_pct - 100.0).abs() < 1e-9);
        kache_free(t);

        assert_eq!(kache_new(ptr::null(), &mut t), KACHE_EINVAL);
        assert!(!kache_strerror(KACHE_EINVAL).is_null());
    }
}
