# Changelog

Notable changes to *k*ache-hash, by release.

## Unreleased

### Added

- `<kache_hash/kache_hash.hpp>`, which includes the whole library.
- `Kmer_Set<k>` and `Kmer_Map<k, T>`, which name the table's concurrent sets and maps.
- `<kache_hash/version.hpp>`, with the version's macros, from which CMake now reads the project's version.
- A compile-time error when AVX2 is missing.
- Installation through CMake: the headers, the libraries, the CMake package `kache_hash` with its components `capi` and `counter`, and the pkg-config files `kache_hash` and `kache_capi`.
- Shared libraries, with `BUILD_SHARED_LIBS`, which export the symbols of their interfaces only.
- `KACHE_HASH_ARCH`, the ISA that the libraries compile for: `native`, the default, `x86-64-v3`, `x86-64-v4`, or `none`.
- `kache_version()` in the C interface, and `version()` in the Rust crate.
- Continuous integration on GitHub Actions: the quick tests with GCC 11–14 and Clang 12–18, in Release and Debug, on Linux and on an Intel Mac; under ASan and UBSan; and for AVX-512 under Intel's SDE; the consumers; and the Rust crates' tests, lints, docs and packages, also on their least Rust, 1.85.
- `TSAN`, a build option for ThreadSanitizer, and a check that the headers', the Rust crates' and `CITATION.cff`'s versions agree.
- Tests for CTest, which `KACHE_HASH_BUILD_TESTS` builds, in top-level builds by default: `quick` ones on small inputs in `test/data/` (248 KB, from `test/data/make-fixtures.py`), `full` ones on the datasets of `scripts/fetch-data.sh small`, and `consumer` ones on installs and sources.
- The Rust crate `kache-hash`, with the safe interface that was `kache-hash-sys`'s, and examples that count the k-mers of a FASTA or FASTQ file from many threads, and tell the k-mers of reads in a reference's strands.
- The crates' metadata for crates.io, and `rust/vendor.sh`, which copies the sources of the C interface's library into `kache-hash-sys`'s package.
- Examples, in `examples/`: in C++, counting k-mers, filtering reads by their k-mers in a reference, strand membership, and counting reads' k-mers within a set with the counter; and in C, strand membership. `KACHE_HASH_BUILD_EXAMPLES` builds them, in top-level builds by default, and CTest checks their outputs on the test inputs against those of `examples/expected/make-expected.py`, a plain count in Python, as CI does the Rust examples'.
- A documentation site, in `docs/`, with guides to installation, the C++, C and Rust interfaces, concurrency, tuning and the supported platforms, and an API reference that Doxygen generates from the headers' comments.
- `REPRODUCE.md`, on reproducing the paper's evaluation of v0.7.0.
- `CONTRIBUTING.md`, with the release checklist, and GitHub's issue templates.

### Changed

- The headers moved to `include/kache_hash/`: include `<kache_hash/Streaming_Kmer_Hash_Table.hpp>`.
  The C interface's and the counter's moved there too: `<kache_hash/kache_capi.h>` and `<kache_hash/Kmer_Set_Counter.hpp>`.
  They need no definitions from the build any more: they set xxHash's themselves, and the cache-line size is 64 bytes unless `KACHE_HASH_CACHE_LINE_SIZE`, or as before `L1_CACHE_LINE_SIZE`, sets it.
- The FASTA and FASTQ parser (`Parser.hpp`, with kseq) and `Directed_Vertex.hpp`, which only the tests and benchmarks use, moved out of the library, to `support/`.
- CMake's options take the prefix `KACHE_HASH_`. The libraries build unless `KACHE_HASH_BUILD_CAPI` or `KACHE_HASH_BUILD_COUNTER` is off, and the benchmarks only with `KACHE_HASH_BUILD_BENCHMARKS`.
- CMake requires GCC 11 or Clang 12 at least, as the headers already did.
- Clang compiles the headers without warnings of GCC's `optimize` attribute.
- `test/check.sh` builds and runs the tests through CMake and CTest.
- `UBSAN` builds stop at undefined behavior, rather than report it and go on.
- The README introduces the library and its installation, and leaves the details to the documentation.
- The crate `kache-hash-sys` holds the raw bindings to the C interface alone; their safe interface moved to the crate `kache-hash`, so `use kache_hash_sys::...` becomes `use kache_hash::...`.
- `kache-hash-sys` compiles the C interface's library with the `cc` crate, without CMake, from its package's copy of the sources or from the repository's tree: for x86-64-v3 by default, for the native CPU if rustc targets it, for x86-64-v4 if rustc's target features include AVX-512's, or for the ISA that `KACHE_HASH_ARCH` names.
  It compiled the library through *k*ache-hash's CMake before, for the native CPU.
- The Rust crates take Rust 1.85 at least.
- *k*ache-hash is now under the BSD 3-Clause License, relicensed from the GNU General Public License v3.0 with the agreement of everyone who has contributed to it.
  Its releases up to v0.7.1 stay under the GPL.

### Deprecated

- The CMake options `COUNTER` and `KACHE_HASH_CAPI`.

### Removed

- The CMake options `AVX_512`, which did nothing, and `CLANG`, which couldn't switch compilers: set `CMAKE_CXX_COMPILER`.
- `scripts/ab/` and `test/check-run.sh`, the A/B and check tooling of the cluster that the paper's runs used.

## 0.7.1 (2026-11-06)

### Changed

- Tables are sized exactly by default, as the paper measures them: the constructor's `pow2_cap`, the C interface's `pow2_capacity`, and the crate's `Options::pow2_capacity` default to false.
  `Kmer_Set_Counter` still rounds its tables up to powers of 2, so that PanGenie's estimate of its k-mers holds them without a resize.

### Added

- `CITATION.cff`.

## 0.7.0 (2026-11-03)

The release that the paper evaluates.

### Added

- A C interface, the library `kache_capi`: concurrent sets of 31-mers, and maps from 31-mers to unsigned values of 8, 16, 32, or 64 bits.
- The Rust crate `kache-hash-sys`, over the C interface.
- Oriented upserts and lookups, which tell a k-mer's occurrences as is from those as its reverse complement: `upsert_oriented`, and `find` with orientations.
- Exact sizing: a table can take the fewest buckets for its k-mers at its load factor, which multiplied high bits map to, rather than a power of 2.
  Tables still round up to powers of 2 by default, and double in rounds.
- Benchmarks: Jellyfish in the paper's experiments, and `bench_kache_hash_exact`, of exactly sized tables.

### Fixed

- Benchmarks: false sharing in the Jellyfish driver, which had slowed it 3–5× at 56 threads.

## 0.6.2 (2026-10-31)

### Changed

- The counter's minimizer size is the tables' default, 22.

### Added

- Benchmarks: a Jellyfish driver, and the measurements of the memory roofs and of the stalls, in `tools/`.

## 0.6.1 (2026-10-31)

### Changed

- The default minimizer size is 22, from 19.

## 0.6.0 (2026-10-30)

### Added

- `update`, which updates the values of the k-mers present in a map, and inserts none, with no tokens or locks.
- `for_each`, a parallel iteration over a table's keys, in shares that partition it.
- Point queries: `find` of one k-mer.
- The counter library, `kache_counter`: `Kmer_Set_Counter`, which counts the k-mers of reads within a set fixed beforehand, as PanGenie counts them with Jellyfish
- `Fastx_Reader`, a parallel reader of FASTA and FASTQ files.
- Benchmarks: the paper's datasets and workloads, with each run's threads and memory placed, and its version and node recorded.

### Changed

- The build configuration works as a subproject's.

## 0.5.1 (2026-10-29)

### Changed

- The overflow table no longer pads its maps' values to their keys' width.

## 0.5.0 (2026-10-29)

### Changed

- Maps pack their slots where the key-value pairs would pad, as with values narrower than keys, so that a value takes only its own bytes.

## 0.4.0 (2026-10-29)

### Changed

- Buckets drop their minimizer coordinates for a spill filter and a lock bit; resizes recompute the minimizers, for faster builds and smaller maps.
- Resizes prefetch their re-inserts in stages, and reserve room for the keys they set aside.
- Inserts count their keys in registers, and check for a pending resize as a run starts.

## 0.3.0 (2026-10-28)

### Changed

- Lookups go through a spill filter per bucket, which marks the minimizers whose keys left it, and through one batch pipeline; their metadata are prefetched and their runs are tracked without branches, for faster queries.
- Without AVX-512, the checksums take a hash of their own.
- Benchmarks: Boost's concurrent tables look k-mers up in bulk, with Boost's default hashing.

## 0.2.0 (2026-10-26)

### Changed

- The streaming operations replace the per-k-mer interface: inserts, upserts, and lookups take whole sequences, or batches of them.
- K-mers are described in batches, with a new minimizer scheme, and on AVX-512 in SIMD lanes.
- The tables sit on huge pages, and the overflow table stores its keys complemented; resizes run in place, and the overflow table grows on its own.
- The streaming operations prefetch ahead, and spin locks pause.

### Added

- Benchmarks: drivers for Boost's tables, libcuckoo, iceberg, CBL and SSHash, their runner and collation, and exact references for their checks.

### Removed

- The example `examples/example.cpp`, and dead code.

## 0.1.1 (2026-10-24)

### Fixed

- GCC 11 miscompiled xxHash's reads at `-O3` with AVX2, breaking the minimizers' hashes: xxHash now reads memory without its aliasing casts.
- A key placed in a bucket's first slot cleared the bucket's lock bit.
- Resizes of sequential tables on several workers lock their buckets, as concurrent tables' do.
- Upserts of k-mers in the overflow table overwrote their values, rather than applied their functions.

### Added

- Benchmarks: a measurement driver with hardware counters, and a script that fetches the datasets.

## 0.1.0 (2026-02-16)

The first release, with the preprint: a streaming k-mer hash table, as sets and maps, with minimizer-based bucketing, concurrent inserts and lookups, and resizes.
