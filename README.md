# *k*ache-hash
*k*ache-hash is a dynamic, concurrent, and cache-efficient hash table for streaming k-mer operations, where sequences of consecutive k-mers are processed in order.
It is highly scalable in throughput with the number of user threads.

## Overview

*k*-mers (substrings of length *k* from DNA sequences) are fundamental building blocks in computational genomics.
Processing billions of *k*-mers efficiently can be accelerated with data structures that can exploit the streaming nature of *k*-mer workloads and modern CPU cache hierarchies.

*k*ache-hash provides a specialized hash table implementation for hash sets and hash maps with *k*-mer keys and arbitrary associated values.
It:

- exploits spatial and temporal locality in streaming k-mer workloads
- Uses minimizer-based bucketing to improve cache efficiency
- supports highly concurrent insertion, lookup, and resize operations
- provides both hash set and hash map interfaces

## Requirements

- C++20 compatible compiler (GCC ≥ 9.1, Clang ≥ 9.0)
- AVX2 support

## Usage

To use *k*ache-hash, copy over the content from its `include` directory to your project, and include the header file `Streaming_Kmer_Hash_Table.hpp` in your code.

The hash table class has the following template signature:

```cpp
template <uint16_t k, bool mt_, typename T_ = void, uint16_t l = 19>
class Streaming_Kmer_Hash_Table
{...}
```

- The `k` parameter defines the *k*-mer size, with `k ≤ 32`.
- The boolean parameter `mt` defines whether the hash table is used in a multi-threaded (`true`) setting or not.
- The `T` parameter defines the type of the associated value, if a *hash map* is used; it is set to `void` otherwise, i.e. a *hash set* is used.
- The `l` parameter defines the minimizer size used in internal hash table operations, with `16 ≤ l < k`.

The table operates on whole sequences: an operation goes through all the *k*-mers of a sequence, or of a batch of sequences such as reads, and characters other than `ACGT` (in either case) split the sequences.
Each user thread registers once for a token, which its updates take:

```cpp
#include "Streaming_Kmer_Hash_Table.hpp"

using namespace kache_hash;

// A concurrent map from 31-mers to their counts, sized for 100M k-mers; it grows as needed.
Streaming_Kmer_Hash_Table<31, true, uint64_t> counts(100'000'000);
const auto token = counts.register_user();  // One per user thread.

counts.upsert(seq, len, [](const uint64_t c){ return c + 1; }, 1, token);   // Count the k-mers of `seq`.
counts.upsert(seqs, lens, n, [](const uint64_t c){ return c + 1; }, 1, token);  // Of `n` sequences at once.
counts.update(seqs, lens, n, [](const uint64_t c){ return c + 1; });   // Only the ones already in the map: no token, no locks.

counts.find(seq, len, [&](const std::size_t i, const std::optional<uint64_t>& c){ /* The k-mer at `seq + i`. */ });
const auto present = counts.count(seq, len);    // Number of k-mers of `seq` in the map.
counts.for_each([&](const Kmer<31>& x, const uint64_t c){ /* Each k-mer of share `w_id` of `w_c`. */ }, w_c, w_id);

// A set of 31-mers.
Streaming_Kmer_Hash_Table<31, true> set(100'000'000);
const auto set_token = set.register_user();
const auto inserted = set.insert(seq, len, set_token);  // Number of k-mers newly inserted.
```

Queries must not overlap updates that can resize the table, and `for_each` must not overlap any updates.
`update`s, which insert nothing, must not overlap other operations but `update`s.

## Benchmarks

The `benchmark` directory has a driver per table: *k*ache-hash, and as baselines, Boost's concurrent and flat tables, libcuckoo, iceberg (on Linux), CBL, and SSHash.
Each builds its table over the *k*-mers of an input and queries it with other inputs, and prints its measurements as TSV.
CBL and SSHash build only with `-DCBL=ON` and `-DSSHASH=ON`, which need `cargo` from `rustup`, and libclang for CBL; SSHash, a static index, is built from GGCAT's maximal unitigs of the input by `scripts/bench/sshash-index.sh`.

```bash
cmake -S . -B build -DCBL=ON -DSSHASH=ON && cmake --build build -j
scripts/fetch-data.sh data small
scripts/bench/run.sh build/benchmark data out chr21 srr5833294-1m ecoli-32
scripts/bench/collate.py out
```

`fetch-data.sh` fetches the small set above, the paper's (`full`, the default; HumGut, no longer served, from a copy, as its header describes), or the earlier campaigns' (`legacy`).
`run.sh` takes its matrix of tables, workloads, sizings, thread counts, and repetitions from the environment, as its header describes.
`collate.py` writes the medians over repetitions to `out/summary.tsv` and the build costs to `out/build-cost.tsv`, and fails if the tables disagree; given several output directories, it collates them as campaigns, e.g. new *k*ache-hash runs against an earlier campaign's baselines.
On explorer, `sbatch scripts/bench/explorer.sbatch <data-dir> <out-dir> <plan>` builds them on a whole node and runs a plan of matrices, such as `scripts/bench/small.plan` or `scripts/bench/full.plan`; `EXPERIMENTS`, `TOOLS`, and `WORKLOADS` select a part of the plan, e.g. `TOOLS=kache-hash` for *k*ache-hash alone.

<!-- ## Citation

If you use this library in academic work, please cite the associated research paper (pending publication). -->
