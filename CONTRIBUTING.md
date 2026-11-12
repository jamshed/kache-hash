# Contributing

Bug reports, questions, and pull requests are welcome.
A bug report helps most with the CPU, the compiler and its flags, *k* and ℓ, the threads, and the commands that reproduce it, which the issue templates ask for.

## Building and testing

```bash
cmake -S . -B build && cmake --build build -j
ctest --test-dir build -L quick -j 8
```

- The tests' labels: `quick`, on the inputs in `test/data/`, in a few minutes; `examples`, the examples on those inputs; `full`, on the datasets of `scripts/fetch-data.sh small` in `data/`; and `consumer`, projects built against installs of *k*ache-hash and against its sources.
- The sanitizers: `-DASAN=ON -DUBSAN=ON`, with Clang, and `-DTSAN=ON`, alone.
- On Apple silicon, `test/check.sh quick` builds for x86-64 and runs the tests under Rosetta.
- The Rust crates, from `rust/`: `cargo test --workspace` and `cargo clippy --workspace --all-targets -- -D warnings`.
- The documentation site: `pip install -r docs/requirements.txt`, Doxygen, and `docs/build.sh`, into `build-docs/site`.

CI runs the quick tests with each supported compiler, the sanitizers, AVX-512 under Intel's SDE, the consumers, and the crates' tests, lints, docs, and packages.

## Style

- Code is indented by 4 spaces, with braces on lines of their own.
- A declaration's documentation is a plain `//` comment before it, of whole sentences, for each overload; Doxygen's filter, `docs/api/doxygen_filter.py`, reads them so. `NB:` marks a note that a reader would otherwise miss.
- Member functions of more than a line are defined out of their class, `inline`.
- Headers include what they use, and nothing else.
- Commit messages are short imperative subjects, such as "Add CTest tests".

## Releasing

1. CI passes on `develop`, and so do the `full` tests, on a CPU with AVX-512 and on one without; an A/B against the previous release, on the paper's indexes, shows no regression beyond noise.
2. The version bumps together in `include/kache_hash/version.hpp`, in `rust/Cargo.toml` (the crates' version, and their dependency on `kache-hash-sys`), and in `CITATION.cff` (`version` and `date-released`): the `version` test checks that they agree.
3. `CHANGELOG.md`'s Unreleased section becomes the release's, with its date.
4. `develop` merges into `main`, and the merge commit takes the tag `vX.Y.Z`, following [SemVer](https://semver.org/).
5. A GitHub Release of the tag carries the changelog's section; Zenodo, where enabled, mints its DOI, which goes into `CITATION.cff` and the README.
6. The crates go to crates.io from the tagged tree, `kache-hash-sys` first, from `rust/`: `./vendor.sh`, `cargo publish -p kache-hash-sys --allow-dirty`, then `cargo publish -p kache-hash --allow-dirty`. `--allow-dirty` admits the vendored copies, which git ignores, so the tree must be clean otherwise.
7. The documentation site deploys from `main`, by the Docs workflow.
8. Bioconda's autobump proposes the release's recipe; its pull request's checks build and test the package.
