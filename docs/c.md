# The C interface

The library `kache_capi` (in `capi/`) puts concurrent sets of 31-mers, and maps from 31-mers to unsigned values of 8, 16, 32, or 64 bits, behind a C interface, for callers in other languages such as Rust.
Its header is `<kache_hash/kache_capi.h>`, and it compiles as C11 and as C++.
A project links it as `kache_hash::kache_capi`, or through pkg-config's `kache_capi`.
CMake builds it unless `KACHE_HASH_BUILD_CAPI` is off.

```c
#include <kache_hash/kache_capi.h>

kache_options opt = kache_default_options(31, 4);   // 31-mers, with 32-bit values.
opt.expected_keys = 100000000;  // Initial sizing. It grows beyond as needed.

kache_table* counts;
kache_new(&opt, &counts);

uint32_t user;
kache_register(counts, &user);  // One registered user per inserting thread.

kache_insert(counts, user, seqs, lens, n, KACHE_ADD, 1, NULL);  // Counts the k-mers of `n` sequences.
kache_count(counts, reads, read_lens, m, present);  // `present[i]`: count of k-mer instances of read `i` in the map.
kache_sum(counts, reads, read_lens, m, sums);   // `sums[i]`: the counts of those k-mers, summed.
kache_lookup(counts, reads, read_lens, m, found, values, NULL); // Each k-mer position's presence and count.
kache_for_each(counts, 1, 0, f, ctx);   // `f(ctx, keys, values, n)` gets the k-mers and their counts, in chunks.
kache_free(counts);

// Forward-strand membership, in a map of 8-bit masks of the orientations each k-mer was inserted in.
kache_insert_oriented(strands, user, seqs, lens, n, KACHE_OR, 1, 2, NULL);
kache_count_masked(strands, reads, read_lens, m, 1, 2, members);    // `members[i]`: read `i`'s k-mers inserted as they occur in it.
```

## Operations

Inserts and updates combine a *k*-mer's value with an operand: keep, assign, add (saturating), or, min, or max.
Oriented inserts choose the operand by whether a *k*-mer occurs as is or as its reverse complement, so that a map of orientation masks gives forward-strand membership.
Queries count, sum, or look up the *k*-mers of batches of sequences, reduce the whole table, or iterate over its keys in shares.
A table's options size it as the C++ constructor does: for `expected_keys` at `load_factor`, exactly, or rounded up to a power of 2 if `pow2_capacity` is 1.
`kache_version()` gives the version of the library, and `KACHE_HASH_VERSION_STRING` that of the header it compiled with.

The functions that can fail return a status code, `KACHE_OK` on success, which `kache_strerror` describes.
Inserts may overlap one another, each user on its own thread, and updates may overlap one another; queries may overlap one another, but neither inserts nor updates.
A table can register at most `KACHE_MAX_USERS` users, 128 for now.

The header documents each function, as does the [API reference](https://jamshed.github.io/kache-hash/api/).
The example [`strand_membership.c`](https://github.com/jamshed/kache-hash/blob/main/examples/strand_membership.c) counts the *k*-mers of reads in a reference's forward and reverse strands.
