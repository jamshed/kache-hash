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
- `kache_version()` in the C interface, and `version()` in the crate.

### Changed

- The headers moved to `include/kache_hash/`: include `<kache_hash/Streaming_Kmer_Hash_Table.hpp>`.
  The C interface's and the counter's moved there too: `<kache_hash/kache_capi.h>` and `<kache_hash/Kmer_Set_Counter.hpp>`.
  They need no definitions from the build any more: they set xxHash's themselves, and the cache-line size is 64 bytes unless `KACHE_HASH_CACHE_LINE_SIZE`, or as before `L1_CACHE_LINE_SIZE`, sets it.
- The FASTA and FASTQ parser (`Parser.hpp`, with kseq) and `Directed_Vertex.hpp`, which only the tests and benchmarks use, moved out of the library, to `support/`.
- CMake's options take the prefix `KACHE_HASH_`. The libraries build unless `KACHE_HASH_BUILD_CAPI` or `KACHE_HASH_BUILD_COUNTER` is off, and the benchmarks only with `KACHE_HASH_BUILD_BENCHMARKS`.
- CMake requires GCC 11 or Clang 12 at least, as the headers already did.
- Clang compiles the headers without warnings of GCC's `optimize` attribute.
- *k*ache-hash is now under the BSD 3-Clause License, relicensed from the GNU General Public License v3.0 with the agreement of everyone who has contributed to it.
  Its releases up to v0.7.1 stay under the GPL.

### Deprecated

- The CMake options `COUNTER` and `KACHE_HASH_CAPI`.

### Removed

- The CMake options `AVX_512`, which did nothing, and `CLANG`, which couldn't switch compilers: set `CMAKE_CXX_COMPILER`.
- `scripts/ab/` and `test/check-run.sh`, the A/B and check tooling of the cluster that the paper's runs used.
