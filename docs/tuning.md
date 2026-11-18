# Tuning

## *k* and the minimizer size ℓ

The tables take *k* ≤ 32, and the C interface, the counter, and the crates *k* = 31.
A *k*-mer's bucket is its ℓ-minimizer's, with 16 ≤ ℓ < *k*, and ℓ = 22 by default: consecutive *k*-mers share their minimizer over runs, and so tend to share a bucket.
In our experience, a larger ℓ speeds up most builds, up to ℓ = 24–26, but slows the negative queries past ℓ = 20.

## Sizing and the load factor

A table sized for `n` *k*-mers takes the fewest buckets of 32 slots that hold them at its load factor, 0.8 by default, and doubles as more are added.
Sizing a table for its distinct *k*-mers, or for an estimate of them, spares its resizes.
In the paper, builds of sets sized for half their *k*-mers, which resize once, ran at 0.65–0.96 of the presized ones' throughput.
A higher load factor saves memory, at the cost of some speed.

Rounding the bucket count up to a power of 2–with the constructor's `pow2_cap`, the C interface's `pow2_capacity`, or `Options::pow2_capacity` in Rust–lowers the load, which speeds the operations, at up to twice the memory of an exact count.
A table has at most `2^32 - 1` buckets, or `2^31` when rounded.

## Threads

Inserts scale with their threads, each with its own token, up to 128 of them.
Queries take no tokens, and support any number of threads.
Resizes run on the table's `resize_worker_c` threads, 16 by default.

## Huge pages

The table maps its buckets, and its overflow table, onto 2 MB pages and advises transparent huge pages for them; this spares TLB misses.
On Linux, `/sys/kernel/mm/transparent_hugepage/enabled` must be `always` or `madvise`; with `never`, the table runs on 4 KB pages.
The benchmarks report the share of their tables on huge pages.

## NUMA

On a machine of several NUMA nodes, interleave the table's memory over them:

```bash
numactl --interleave=all ./app
```

An interleave takes free memory on every node: a node full of page cache, such as one that read the input, takes fewer of the table's pages, and the table runs slower on the skew.
Dropping the page cache before a run, or leaving room on every node, keeps the interleave even.

## The allocator

The table maps its buckets itself, so the allocator matters little to it.
The paper's benchmarks link jemalloc 5.3.0 into every table, for the baselines' allocations.

## The ISA

The compiled libraries target the building CPU by default, and the header-only table the consumer's flags.
