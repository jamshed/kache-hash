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
