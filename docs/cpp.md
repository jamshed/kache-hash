# The C++ tables

The table is the class template `kache_hash::Streaming_Kmer_Hash_Table`, in `<kache_hash/Streaming_Kmer_Hash_Table.hpp>`:

```cpp
template <uint16_t k, bool mt_, typename T_ = void, uint16_t l = 22>
class Streaming_Kmer_Hash_Table
{...}
```

- The `k` parameter defines the *k*-mer size, with `k ≤ 32`.
- The boolean parameter `mt` defines whether the hash table is used in a multi-threaded (`true`) setting or not.
- The `T` parameter defines the type of the associated value, if a *hash map* is used; it is set to `void` otherwise, i.e. a *hash set* is used.
- The `l` parameter defines the minimizer size used in internal hash table operations, with `16 ≤ l < k`.

`Kmer_Set<k>` and `Kmer_Map<k, T>`, in `<kache_hash/kache_hash.hpp>`, name its concurrent sets and maps, with `l = 22` unless a last parameter gives it.

## Operations

The table operates on whole sequences: an operation goes through all the *k*-mers of a sequence, or of a batch of sequences such as reads.
Characters other than `ACGT` split the sequences.
A table keeps each *k*-mer in its canonical form, the lesser of it and its reverse complement, so that lookups find a *k*-mer in either orientation.
Each user thread registers once for a *token*, which its inserts require:

```cpp
#include <kache_hash/Streaming_Kmer_Hash_Table.hpp>

using namespace kache_hash;

// A concurrent map from 31-mers to their counts, sized for 100M k-mers. It grows as needed.
Streaming_Kmer_Hash_Table<31, true, uint64_t> counts(100'000'000);
const auto token = counts.register_user();  // One per user thread.

counts.upsert(seq, len, [](const uint64_t c){ return c + 1; }, 1, token);   // Count the k-mers of `seq`.
counts.upsert(seqs, lens, n, [](const uint64_t c){ return c + 1; }, 1, token);  // Of `n` sequences at once.
counts.update(seqs, lens, n, [](const uint64_t c){ return c + 1; });   // Only the ones already in the map: no token, no locks.

counts.find(seq, len, [&](const std::size_t i, const std::optional<uint64_t>& c){ /* The k-mer at `seq + i`. */ });
const std::optional<uint64_t> c = counts.find(kmer);  // One k-mer: the 31 characters at `kmer`.
const auto present = counts.count(seq, len);    // Number of k-mers of `seq` in the map.
counts.for_each([&](const Kmer<31>& x, const uint64_t c){ /* Each k-mer of worker `w_id`'s share, out of `w_c`. */ }, w_c, w_id);

// A set of 31-mers.
Streaming_Kmer_Hash_Table<31, true> set(100'000'000);
const auto set_token = set.register_user();
const auto inserted = set.insert(seq, len, set_token);  // Number of k-mers newly inserted.
```

Batching sequences, such as reads, into one call lets the table describe their *k*-mers at once.
This executes faster than a call per sequence.
A map's `upsert` gives an absent *k*-mer the value it takes, and applies its function to a present one's.
A map's `update` applies the function only to the present *k*-mers, and inserts none.
`for_each` visits the *k*-mers in the share of worker `w_id`, of `w_c` shares; the shares partition the table, so that `w_c` threads can visit each *k*-mer once between them.

## Orientations

`upsert_oriented` and the batch `find` can tell a *k*-mer's occurrences as is, i.e. as its canonical form, from those as its reverse complement:

```cpp
// A map from 31-mers to the orientations they occur in: bit 1 as their canonical form, bit 2 as its reverse complement.
Streaming_Kmer_Hash_Table<31, true, uint8_t> strands(100'000'000);
const auto strand_token = strands.register_user();
strands.upsert_oriented(seqs, lens, n, [](const uint8_t o, const bool fwd){ return uint8_t(o | (fwd ? 1 : 2)); }, strand_token);   // `fwd`: whether it occurs as its canonical form.
strands.find(seqs, lens, n, [&](const std::size_t q, const std::size_t i, const std::optional<uint8_t>& o, const bool fwd){ /* The k-mer at `seqs[q] + i`. */ });
```

A read's *k*-mer then occurs in the forward strand of the sequences inserted if it occurs there in its orientation in the read, as the example [`strand_membership.cpp`](https://github.com/jamshed/kache-hash/blob/main/examples/strand_membership.cpp) counts.

## Sizing

A table sized for `n` *k*-mers allocates the fewest buckets of 32 entries that hold them at its load factor, 0.8 by default, and doubles as more are added.
After `n`, the constructor takes the number of threads for resizes, the load factor, and whether to round the bucket count up to a power of 2, false by default.
A power of 2 sizing lowers the load, which speeds the operations, at up to twice the memory at most.
A table has at most `2^32 - 1` buckets; `fits` tells whether a size is within that.
`size`, `capacity`, `bucket_bytes`, and `key_distribution` describe a table: its *k*-mers, its slots, the bytes of its buckets, and the shares of its keys in their primary buckets, their secondary buckets, and the overflow table.

## Concurrency

Queries must not overlap inserts that can grow the table, and `for_each` must not overlap a modification.
`update`s, which insert nothing, must not overlap other operations but `update`s.
[Concurrency](concurrency.md) tabulates which operations may overlap.

The [API reference](https://jamshed.github.io/kache-hash/api/) documents each operation.
