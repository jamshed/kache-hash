# kache-hash

Concurrent tables of *k*-mers from [*k*ache-hash](https://github.com/jamshed/kache-hash), a dynamic, concurrent and cache-efficient hash table for streaming *k*-mer operations: sets of 31-mers, and maps from 31-mers to `u8`, `u16`, `u32` or `u64` values, built and queried from many threads.

A table is built through a `Builder`, with an `Inserter` per thread, then frozen into a `Table` for queries from any number of threads; a map's `Updater` borrows it to update it.
So the borrow checker keeps inserts, updates and queries from overlapping, as the table requires.
Inserts and updates combine a k-mer's value with an operand: keep, assign, add (saturating), or, min, or max.
Queries count, sum, or look up the k-mers of batches of sequences, reduce the whole table, or iterate over its keys in shares.
A table keeps each k-mer in its canonical form, so that it finds k-mers in either orientation; oriented inserts and masked counts tell a k-mer's occurrences as is from those as its reverse complement, e.g. for strand membership.

```rust
use kache_hash::{Builder, Op, Options};

let b = Builder::<u32>::new(&Options::new(31))?;    // A map from 31-mers to 32-bit counts.
b.inserter()?.insert(&seqs, Op::Add, 1)?;   // An inserter per thread.
let counts = b.freeze();    // Queries from any number of threads.
counts.count(&reads, &mut present)?;    // `present[i]`: the k-mer positions of `reads[i]` in the map.
```

The examples count the k-mers of a FASTA or FASTQ file from a number of threads, and tell the k-mers of reads in a reference's forward and reverse strands:

```bash
cargo run --release --example kmer_counts -- seqs.fa 8
cargo run --release --example strand_membership -- reference.fa reads.fq
```

## Building

The crate [`kache-hash-sys`](https://crates.io/crates/kache-hash-sys) compiles *k*ache-hash's C interface library with the `cc` crate: it takes a C++20 compiler, GCC 11 or Clang 12, or later, and an x86-64 CPU with AVX2.
A default build runs on any such CPU; `RUSTFLAGS="-C target-cpu=native"` builds the fastest tables for the CPU that builds them, and `KACHE_HASH_ARCH` picks an ISA, as `kache-hash-sys` documents.

## Citation

If you use *k*ache-hash, please cite its paper:

> J. Khan, R. Patro, and P. Pandey. *k*ache-hash: A Dynamic, Concurrent, and Cache-Efficient Hash Table for Streaming *k*-mer Operations. bioRxiv, 2026. [doi:10.64898/2026.02.13.705625](https://doi.org/10.64898/2026.02.13.705625)

## License

BSD-3-Clause, as *k*ache-hash's.
