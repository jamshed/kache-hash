// Checks the crate against exact references: forward-strand membership, each operation's semantics, concurrent
// counting from a tiny table that resizes, a set's queries across a thaw, concurrent updates, sizing, and errors and
// panics.
// The sequences are random, with an `N` every 997th base, some reverse-complemented.

use kache_hash::{supports, Builder, Error, Op, Options, Reduction, ABSENT, NOT_KMER, PRESENT};
use std::collections::{HashMap, HashSet};

const K: usize = 31;

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            x => x,
        })
        .collect()
}

/// A deterministic random DNA sequence of `n` bases, with an `N` every 997th.
fn random_seq(n: usize, seed: &mut u64) -> Vec<u8> {
    (0..n)
        .map(|i| {
            *seed ^= *seed << 13;
            *seed ^= *seed >> 7;
            *seed ^= *seed << 17;
            if i % 997 == 996 {
                b'N'
            } else {
                b"ACGT"[(*seed >> 32) as usize & 3]
            }
        })
        .collect()
}

/// The 2-bit encoding of the k-mer `w`, with A, C, G and T as 0 to 3 and the first base the most significant; `None`
/// if it has an `N`.
fn encode(w: &[u8]) -> Option<u64> {
    w.iter().try_fold(0u64, |x, &c| match c {
        b'A' => Some(x << 2),
        b'C' => Some(x << 2 | 1),
        b'G' => Some(x << 2 | 2),
        b'T' => Some(x << 2 | 3),
        _ => None,
    })
}

/// The canonical encoding of each k-mer position of `s`, and whether it occurs as is; `None` for those with an `N`.
fn kmers(s: &[u8]) -> Vec<Option<(u64, bool)>> {
    s.windows(K)
        .map(|w| {
            let f = encode(w)?;
            let r = encode(&revcomp(w)).unwrap();
            Some((f.min(r), f <= r))
        })
        .collect()
}

fn sequences(c: usize, seed: &mut u64) -> Vec<Vec<u8>> {
    (0..c)
        .map(|i| {
            let s = random_seq(300 + (i * 37) % 1500, seed);
            if i % 3 == 0 {
                revcomp(&s)
            } else {
                s
            }
        })
        .collect()
}

fn options(expected: u64, workers: usize) -> Options {
    Options { expected_keys: expected, resize_workers: workers, ..Options::new(K) }
}

#[test]
fn support() {
    assert!(supports::<()>(31) && supports::<u8>(31) && supports::<u64>(31));
    assert!(!supports::<u8>(25) && !supports::<()>(32));
}

#[test]
fn forward_strand_membership() {
    let x = b"ACGTTGCATGTCGCATGATGCATGAGAGCTA";
    let r = revcomp(x);

    let b = Builder::<u8>::new(&options(1000, 1)).unwrap();
    b.inserter().unwrap().insert_oriented(&[x], Op::Or, 1, 2).unwrap();
    let m = b.freeze();
    let mut c = [0u64; 2];
    m.count_masked(&[&x[..], &r[..]], 1, 2, &mut c).unwrap();
    assert_eq!(c, [1, 0]);
    m.count(&[&x[..], &r[..]], &mut c).unwrap();
    assert_eq!(c, [1, 1]);
    assert_eq!(m.reduce(Reduction::Popcount, 3, 2).unwrap(), 1);

    let b = m.thaw();
    b.inserter().unwrap().insert_oriented(&[&r[..]], Op::Or, 1, 2).unwrap();
    let m = b.freeze();
    m.count_masked(&[&x[..], &r[..]], 1, 2, &mut c).unwrap();
    assert_eq!(c, [1, 1]);
    assert_eq!(m.reduce(Reduction::Popcount, 3, 2).unwrap(), 2);
    assert_eq!((m.len(), m.stats().unwrap().keys), (1, 1));
}

#[test]
fn operations() {
    let x: &[u8] = b"ACGTTGCATGTCGCATGATGCATGAGAGCTA";
    for (op, a, b, want) in [(Op::Keep, 5u8, 9, 5), (Op::Assign, 5, 9, 9), (Op::Add, 200, 100, 255), (Op::Or, 1, 4, 5), (Op::Min, 5, 9, 5), (Op::Max, 5, 9, 9)] {
        let bl = Builder::<u8>::new(&options(1000, 1)).unwrap();
        let mut ins = bl.inserter().unwrap();
        assert_eq!(ins.insert(&[x], op, a).unwrap(), 1);
        assert_eq!(ins.insert(&[x], op, b).unwrap(), 0);
        drop(ins);
        let m = bl.freeze();
        assert_eq!(m.get(x).unwrap(), Some(want), "{op:?}");
        assert_eq!(m.get(&revcomp(x)).unwrap(), Some(want), "{op:?}, reverse complement");
        assert_eq!(m.get(b"ACGTTGCATGTCGCATGATGCATGAGAGCTN").unwrap(), None);
    }

    // Oriented inserts take the operand by orientation, and refuse the operations an absent k-mer can't start from.
    let bl = Builder::<u16>::new(&options(1000, 1)).unwrap();
    let mut ins = bl.inserter().unwrap();
    ins.insert_oriented(&[x], Op::Add, 1, 1000).unwrap();
    ins.insert_oriented(&[&revcomp(x)[..]], Op::Add, 1, 1000).unwrap();
    assert_eq!(ins.insert_oriented(&[x], Op::Keep, 1, 2), Err(Error::INVALID));
    assert_eq!(ins.insert_oriented(&[x], Op::Min, 1, 2), Err(Error::INVALID));
    drop(ins);
    assert_eq!(bl.freeze().get(x).unwrap(), Some(1001));
}

#[test]
fn concurrent_counts_against_a_reference() {
    let mut seed = 0x9E37_79B9_7F4A_7C15u64;
    let seqs = sequences(400, &mut seed);
    let mut reference: HashMap<u64, u32> = HashMap::new();
    for s in &seqs {
        for (x, _) in kmers(s).into_iter().flatten() {
            *reference.entry(x).or_default() += 1;
        }
    }

    let b = Builder::<u32>::new(&options(10, 4)).unwrap();
    std::thread::scope(|sc| {
        for t in 0..4 {
            let (b, seqs) = (&b, &seqs);
            sc.spawn(move || {
                let mut ins = b.inserter().unwrap();
                let mine: Vec<&[u8]> = seqs.iter().skip(t).step_by(4).map(|s| s.as_slice()).collect();
                for chunk in mine.chunks(7) {
                    ins.insert(chunk, Op::Add, 1).unwrap();
                }
            });
        }
    });
    let m = b.freeze();
    assert_eq!(m.len(), reference.len() as u64);
    assert_eq!(m.reduce(Reduction::Sum, u32::MAX, 3).unwrap(), reference.values().map(|&v| v as u64).sum::<u64>());

    // Every key once, with its count, over 3 shares iterated concurrently.
    let shares: Vec<Vec<(u64, u32)>> = std::thread::scope(|sc| {
        let hs: Vec<_> = (0..3)
            .map(|s| {
                let m = &m;
                sc.spawn(move || {
                    let mut kv = Vec::new();
                    m.for_each_chunk(3, s, |keys, vals| {
                        kv.extend(keys.iter().copied().zip(vals.iter().copied()));
                        true
                    })
                    .unwrap();
                    kv
                })
            })
            .collect();
        hs.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let seen: HashMap<u64, u32> = shares.iter().flatten().copied().collect();
    assert_eq!(shares.iter().map(|s| s.len()).sum::<usize>(), reference.len());
    assert_eq!(seen, reference);

    // The inserted sequences, their reverse complements, and new ones: per-sequence counts and sums, and per-position
    // lookups.
    let mut queries = seqs.clone();
    queries.extend(seqs.iter().step_by(5).map(|s| revcomp(s)));
    queries.extend((0..50).map(|_| random_seq(2000, &mut seed)));
    queries.extend([Vec::new(), b"ACGT".to_vec(), vec![b'N'; 40]]);
    let (mut counts, mut sums) = (vec![0u64; queries.len()], vec![0u64; queries.len()]);
    m.count(&queries, &mut counts).unwrap();
    m.sum(&queries, &mut sums).unwrap();
    let p = m.positions(&queries);
    let (mut found, mut vals, mut fwd) = (vec![9u8; p], vec![7u32; p], vec![9u8; p]);
    m.lookup_values(&queries, &mut found, &mut vals, Some(&mut fwd)).unwrap();

    let mut o = 0;
    for (i, q) in queries.iter().enumerate() {
        let (mut c, mut s) = (0, 0);
        for w in kmers(q) {
            match w {
                None => assert_eq!((found[o], vals[o]), (NOT_KMER, 0)),
                Some((x, f)) => {
                    let v = reference.get(&x).copied();
                    assert_eq!((found[o], vals[o], fwd[o] != 0), (if v.is_some() { PRESENT } else { ABSENT }, v.unwrap_or(0), f));
                    c += v.is_some() as u64;
                    s += v.unwrap_or(0) as u64;
                }
            }
            o += 1;
        }
        assert_eq!((counts[i], sums[i]), (c, s), "query {i}");
    }
    assert_eq!(o, p);
}

#[test]
fn a_set_across_a_thaw() {
    let mut seed = 0x0123_4567_89AB_CDEFu64;
    let seqs = sequences(60, &mut seed);
    let (first, second) = seqs.split_at(30);
    let reference: HashSet<u64> = first.iter().flat_map(|s| kmers(s)).flatten().map(|(x, _)| x).collect();

    let b = Builder::<()>::new(&options(10, 2)).unwrap();
    std::thread::scope(|sc| {
        for half in first.chunks(15) {
            let b = &b;
            sc.spawn(move || b.inserter().unwrap().insert(half).unwrap());
        }
    });
    let s = b.freeze();
    assert_eq!(s.len(), reference.len() as u64);
    let mut keys = HashSet::new();
    s.for_each_chunk(1, 0, |k, v| {
        assert_eq!(k.len(), v.len());
        keys.extend(k.iter().copied());
        true
    })
    .unwrap();
    assert_eq!(keys, reference);

    let q = [&second[0][..], &first[0][..]];
    let mut found = vec![0u8; s.positions(&q)];
    s.lookup(&q, &mut found, None).unwrap();
    let want: Vec<u8> = q.iter().flat_map(|s| kmers(s)).map(|w| w.map_or(NOT_KMER, |(x, _)| if reference.contains(&x) { PRESENT } else { ABSENT })).collect();
    assert_eq!(found, want);
    assert!(s.contains(&first[0][..K]).unwrap() || first[0][..K].contains(&b'N'));

    // More inserts, through a released user.
    let b = s.thaw();
    let inserted = b.inserter().unwrap().insert(second).unwrap();
    let s = b.freeze();
    let all: HashSet<u64> = seqs.iter().flat_map(|s| kmers(s)).flatten().map(|(x, _)| x).collect();
    assert_eq!((s.len(), inserted), (all.len() as u64, (all.len() - reference.len()) as u64));
}

#[test]
fn concurrent_updates() {
    let mut seed = 0xDEAD_BEEF_0BAD_F00Du64;
    let seqs = sequences(100, &mut seed);
    let (index, queries) = seqs.split_at(50);
    let b = Builder::<u16>::new(&options(1 << 16, 2)).unwrap();
    b.inserter().unwrap().insert(index, Op::Assign, 0).unwrap();
    let mut m = b.freeze();
    let keys: HashSet<u64> = index.iter().flat_map(|s| kmers(s)).flatten().map(|(x, _)| x).collect();

    // Each present k-mer of the queries and the index gains 1 per occurrence; absent ones stay absent.
    let all: Vec<&[u8]> = queries.iter().chain(index).map(|s| s.as_slice()).collect();
    let up = m.updater();
    let present: u64 = std::thread::scope(|sc| {
        let up = &up;
        let hs: Vec<_> = all.chunks(13).map(|c| sc.spawn(move || up.update(c, Op::Add, 1).unwrap())).collect();
        hs.into_iter().map(|h| h.join().unwrap()).sum()
    });
    let mut reference: HashMap<u64, u16> = HashMap::new();
    for s in &all {
        for (x, _) in kmers(s).into_iter().flatten().filter(|(x, _)| keys.contains(x)) {
            *reference.entry(x).or_default() += 1;
        }
    }
    assert_eq!(present, reference.values().map(|&v| v as u64).sum::<u64>());
    assert_eq!(m.len(), keys.len() as u64);
    m.for_each_chunk(1, 0, |k, v| {
        for (x, v) in k.iter().zip(v) {
            assert_eq!(*v, reference.get(x).copied().unwrap_or(0));
        }
        true
    })
    .unwrap();
}

#[test]
fn sizing() {
    // Exactly for the expected k-mers at the load factor, or rounded up to a power of 2.
    for (pow2_capacity, capacity) in [(false, 1280), (true, 2048)] {
        let b = Builder::<()>::new(&Options { pow2_capacity, ..options(1000, 1) }).unwrap();
        assert_eq!(b.freeze().stats().unwrap().capacity, capacity);
    }

    // Beyond the most buckets.
    assert_eq!(Builder::<()>::new(&Options { expected_keys: 1 << 37, pow2_capacity: false, ..Options::new(K) }).err(), Some(Error::INVALID));
    assert_eq!(Builder::<()>::new(&Options { expected_keys: 55 << 30, pow2_capacity: true, ..Options::new(K) }).err(), Some(Error::INVALID));
}

#[test]
fn errors_and_panics() {
    assert_eq!(Builder::<u8>::new(&Options { k: 25, ..Options::new(K) }).err(), Some(Error::INVALID));
    assert_eq!(Builder::<u8>::new(&Options { load_factor: 1.0, ..Options::new(K) }).err(), Some(Error::INVALID));
    assert!(Error::INVALID.to_string().contains("invalid argument"));

    let mut seed = 7u64;
    let seqs = sequences(200, &mut seed);
    let b = Builder::<u64>::new(&options(1 << 16, 2)).unwrap();
    b.inserter().unwrap().insert(&seqs, Op::Add, 1).unwrap();
    let m = b.freeze();

    // An iteration that stops gets one chunk; one that panics passes the panic on.
    let mut chunks = 0;
    m.for_each_chunk(1, 0, |_, _| {
        chunks += 1;
        false
    })
    .unwrap();
    assert_eq!(chunks, 1);
    let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| m.for_each_chunk(1, 0, |_, _| panic!("in a chunk"))));
    assert!(r.is_err());
    assert_eq!(m.for_each_chunk(2, 2, |_, _| true), Err(Error::INVALID));
    assert_eq!(m.reduce(Reduction::Keys, u64::MAX, 0), Err(Error::INVALID));
}

#[test]
fn version() {
    // The crate and its C library release together.
    assert_eq!(kache_hash::version(), env!("CARGO_PKG_VERSION"));
}
