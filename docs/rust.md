# The Rust crates

Two crates, in `rust/`, bring the C interface to Rust:

- [`kache-hash`](https://docs.rs/kache-hash), the safe interface: typed sets and maps of `u8` to `u64` values;
- [`kache-hash-sys`](https://docs.rs/kache-hash-sys), the raw bindings to `kache_capi.h`, and the build of its library.

```bash
cargo add kache-hash
```

A table is built through a `Builder`, with an `Inserter` per thread, then frozen into a `Table` for queries from any number of threads; a map's `Updater` borrows it to update it, so that the borrow checker keeps the phases apart.

```rust
use kache_hash::{Builder, Op, Options};

let b = Builder::<u32>::new(&Options::new(31))?;    // A map from 31-mers to 32-bit counts.
b.inserter()?.insert(&seqs, Op::Add, 1)?;   // An inserter per thread.
let counts = b.freeze();    // Queries from any number of threads.
counts.count(&reads, &mut present)?;    // `present[i]`: the k-mer positions of `reads[i]` in the map.
```

## The build

`kache-hash-sys` compiles the C interface's library with the [`cc`](https://crates.io/crates/cc) crate, with no CMake: from the sources in its package, or from *k*ache-hash's tree when the crate builds within it.
It takes GCC 11 or Clang 12, or later, as `CXX` or the system's `c++`, and Rust 1.85 or later.
The library compiles with optimizations in every Cargo profile, for:

- the ISA that `KACHE_HASH_ARCH` names, if set: `native`, `x86-64-v3`, `x86-64-v4`, or `none`, for the compiler's own flags, such as those of `CXXFLAGS`;
- else, the native CPU, if rustc targets it, e.g. with `-C target-cpu=native`;
- else, x86-64-v4 (AVX-512), if rustc's target features include AVX-512's F, BW, CD, DQ, and VL, e.g. with `-C target-cpu=x86-64-v4`;
- and x86-64-v3 (AVX2) otherwise.

So a default build runs on any x86-64 CPU with AVX2, and the fastest tables for the CPU that builds them come from

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

## Examples

The crate's examples count the *k*-mers of a FASTA or FASTQ file from a number of threads, and tell the *k*-mers of reads in a reference's forward and reverse strands:

```bash
cd rust
cargo run --release --example kmer_counts -- seqs.fa 8
cargo run --release --example strand_membership -- reference.fa reads.fq
```
