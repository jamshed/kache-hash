# Benchmarks

The `benchmark` directory has a driver per table: *k*ache-hash, and as baselines, Boost's concurrent and flat tables, libcuckoo, iceberg (on Linux), Jellyfish, CBL, and SSHash.
*k*ache-hash's driver rounds its tables' bucket counts up to powers of 2, and a second one, `bench_kache_hash_exact` (the table `kache-hash-exact` to `run.sh`), sizes them exactly, as the tables do by default.
Each builds its table over the *k*-mers of an input and queries it with other inputs, and prints its measurements as TSV.
CMake builds them with `-DKACHE_HASH_BUILD_BENCHMARKS=ON`, fetching the baselines; CBL and SSHash only with `-DCBL=ON` and `-DSSHASH=ON` too, which need `cargo` from `rustup`, and libclang for CBL; SSHash, a static index, is built from GGCAT's maximal unitigs of the input by `scripts/bench/sshash-index.sh`.

```bash
cmake -S . -B build -DKACHE_HASH_BUILD_BENCHMARKS=ON -DCBL=ON -DSSHASH=ON && cmake --build build -j
scripts/fetch-data.sh data small
scripts/bench/run.sh build/benchmark data out chr21 srr5833294-1m ecoli-32
scripts/bench/collate.py out
```

`fetch-data.sh` fetches the small set above, the paper's (`full`, the default; HumGut, no longer served, from a copy, as its header describes), or the earlier campaigns' (`legacy`).
`run.sh` takes its matrix of tables, workloads, sizings, thread counts, and repetitions from the environment, as its header describes.
`collate.py` writes the medians over repetitions to `out/summary.tsv` and the build costs to `out/build-cost.tsv`, and fails if the tables disagree; given several output directories, it collates them as campaigns, e.g. new *k*ache-hash runs against an earlier campaign's baselines.
On explorer, `sbatch scripts/bench/explorer.sbatch <data-dir> <out-dir> <plan>` builds them on a whole node and runs a plan of matrices, such as `scripts/bench/small.plan` or `scripts/bench/full.plan`; `EXPERIMENTS`, `TOOLS`, and `WORKLOADS` select a part of the plan, e.g. `TOOLS=kache-hash` for *k*ache-hash alone.
Its `#SBATCH` lines and environment are explorer's, as the paper's runs took them; its header tells what to change for another Slurm cluster.

## The paper's evaluation

The paper evaluates *k*ache-hash v0.7.0 with these benchmarks, on whole nodes of two Intel Xeon Platinum 8276 (Cascade Lake).
`scripts/bench/paper.plan` lists its experiments, each a `run.sh` matrix under its name, with comments on what it measures; `explorer.sbatch` and `scripts/bench/feed.sh` ran them, and `collate.py` collated their outputs.
