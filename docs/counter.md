# Counting *k*-mers within a set

The library `kache_counter`, in `counter/`, counts the *k*-mers of reads that belong to a set fixed beforehand.
The *k*-mer set is frozen here to match [PanGenie](https://github.com/eblerjana/pangenie)'s *k*-mer counting workload with Jellyfish, on a frozen table of 32-bit counts.
Its interface is plain C++, without the table's headers.
Therefore a project linking it as `kache_hash::kache_counter` keeps its own compile flags.
CMake builds it unless `KACHE_HASH_BUILD_COUNTER` is off.

```cpp
#include <kache_hash/Kmer_Set_Counter.hpp>

kache_hash::Kmer_Set_Counter counter(31, 3'000'000'000, 24);   // 31-mers; room for 3G of them; 24 user threads.
counter.prime("segments.fa");   // The set: the k-mers of these sequences, set at count 0.
const auto found = counter.count("reads.fq");   // Counts the reads' k-mers in the set; returns the occurrences counted.
const auto bins = counter.histogram(10'000);    // `bins[c]`: the number of k-mers with count `c`.
const auto c = counter.count_of(kmer);  // One k-mer's count: the 31 characters at `kmer`.
```

The counter takes 31-mers only, *for now*.
Its table rounds its bucket count up to a power of 2, at load factor 0.9 by default.
Its operations must not overlap one another, except for `histogram`s and `count_of`s.

## Reading FASTA and FASTQ in parallel

The library's `Fastx_Reader`, in `<kache_hash/Fastx_Reader.hpp>`, reads the sequences of a FASTA or FASTQ file in parallel, and feeds them to the counter's threads.
The first character sets the format: `>` for FASTA, `@` for FASTQ.
A regular file is memory-mapped and split among the workers; any other input, such as a pipe, is read by one more thread, which hands them chunks of whole records.

```cpp
#include <kache_hash/Fastx_Reader.hpp>

kache_hash::Fastx_Reader reader("reads.fq");
reader.read(16, [&](const std::size_t t, const char* const* seqs, const std::size_t* lens, const std::size_t c)
{
    // Worker `t` gets a batch of `c` sequences, which a table's batch operations take as they come.
});
```

The example [`count_in_set.cpp`](https://github.com/jamshed/kache-hash/blob/main/examples/count_in_set.cpp) counts the *k*-mers of reads in a reference's.
