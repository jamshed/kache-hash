# kache-hash's API reference

The reference of kache-hash's C++ and C interfaces, from the comments of their headers; [the guides](https://jamshed.github.io/kache-hash/) introduce them, and [docs.rs](https://docs.rs/kache-hash) documents the Rust crates.

- **The table**: the class template `kache_hash::Streaming_Kmer_Hash_Table`, in `<kache_hash/Streaming_Kmer_Hash_Table.hpp>`, and its concurrent sets and maps, `kache_hash::Kmer_Set` and `kache_hash::Kmer_Map`; `<kache_hash/kache_hash.hpp>` includes them all.
- **The counter**: `kache_hash::Kmer_Set_Counter`, which counts the k-mers of reads within a set fixed beforehand, and its parallel reader of FASTA and FASTQ files, `kache_hash::Fastx_Reader`, in the library `kache_counter`.
- **The C interface**: `kache_capi.h`, the functions `kache_*` of the library `kache_capi`.
- **Configuration**: `config.hpp`, what the headers take from the compiler, and `version.hpp`, the version's macros.
