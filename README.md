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

<!-- ## Citation

If you use this library in academic work, please cite the associated research paper (pending publication). -->
