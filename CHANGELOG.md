# Changelog

Notable changes to *k*ache-hash, by release.

## Unreleased

### Added

- `<kache_hash/kache_hash.hpp>`, which includes the whole library.
- `Kmer_Set<k>` and `Kmer_Map<k, T>`, which name the table's concurrent sets and maps.
- `<kache_hash/version.hpp>`, with the version's macros, from which CMake now reads the project's version.
- A compile-time error when AVX2 is missing.

### Changed

- The headers moved to `include/kache_hash/`: include `<kache_hash/Streaming_Kmer_Hash_Table.hpp>`.
  They need no definitions from the build any more: they set xxHash's themselves, and the cache-line size is 64 bytes unless `KACHE_HASH_CACHE_LINE_SIZE`, or as before `L1_CACHE_LINE_SIZE`, sets it.
- The FASTA and FASTQ parser (`Parser.hpp`, with kseq) and `Directed_Vertex.hpp`, which only the tests and benchmarks use, moved out of the library, to `support/`.
- *k*ache-hash is now under the BSD 3-Clause License, relicensed from the GNU General Public License v3.0 with the agreement of everyone who has contributed to it.
  Its releases up to v0.7.1 stay under the GPL.
