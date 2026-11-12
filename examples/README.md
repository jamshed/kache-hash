# Examples

Small programs on *k*ache-hash's tables of 31-mers, in C++, C, and Rust.
They read uncompressed FASTA or FASTQ files and print plain text, so their outputs can be compared across the languages.

| Example | Language | What it does | Demonstrates |
|---|---|---|---|
| [`count_kmers.cpp`](count_kmers.cpp) | C++ | Counts the k-mers of a file from a number of threads, and prints their number, their occurrences, and a histogram of their counts | `Kmer_Map`, a token per thread, `upsert`, `for_each` in shares |
| [`filter_reads.cpp`](filter_reads.cpp) | C++ | Keeps the reads with at least a share of their k-mers in a reference, as [K2Rmini](https://github.com/jamshed/K2Rmini) does, in either strand | `Kmer_Set`, `insert`, concurrent `count`s |
| [`strand_membership.cpp`](strand_membership.cpp) | C++ | Counts each read's k-mers in a reference's forward and reverse strands | `upsert_oriented`, `find` with orientations |
| [`count_in_set.cpp`](count_in_set.cpp) | C++ | Counts the occurrences in reads of a reference's k-mers, as [PanGenie](https://github.com/jamshed/pangenie) counts its graph's | The counter library, `Kmer_Set_Counter` |
| [`strand_membership.c`](strand_membership.c) | C | As `strand_membership.cpp`, through the C interface | `kache_insert_oriented`, `kache_count_masked` |
| [`kmer_counts.rs`](../rust/kache-hash/examples/kmer_counts.rs) | Rust | As `count_kmers.cpp` | `Builder`, an `Inserter` per thread, `for_each_chunk` in shares |
| [`strand_membership.rs`](../rust/kache-hash/examples/strand_membership.rs) | Rust | As `strand_membership.cpp` | `insert_oriented`, `count_masked` |

The C++ and C examples build with *k*ache-hash's top-level CMake build, unless `-DKACHE_HASH_BUILD_EXAMPLES=OFF`, into `examples/` of the build directory.
[`fastx.hpp`](fastx.hpp) is the C++ examples' reader of input files.
The Rust ones build and run from `rust/`:

```bash
cmake -S . -B build && cmake --build build -j
gzip -dc test/data/ecoli.fa.gz > ecoli.fa && gzip -dc test/data/reads.fq.gz > reads.fq
build/examples/count_kmers ecoli.fa 8
build/examples/filter_reads ecoli.fa reads.fq 0.5 8 > kept.fa
(cd rust && cargo run --release --example strand_membership -- ../ecoli.fa ../reads.fq)
```

CTest runs the C++ and C examples on the inputs in `test/data/`, as the tests labeled `examples`, and compares their outputs with those in [`expected/`](expected), which [`make-expected.py`](expected/make-expected.py) makes by a plain count of the inputs' canonical 31-mers.
CI runs the Rust examples against the same outputs.
