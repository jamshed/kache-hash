# kache-hash-sys

Raw bindings to the C interface of [*k*ache-hash](https://github.com/jamshed/kache-hash), a concurrent and cache-efficient hash table for streaming *k*-mer operations, and the build of its library.
The crate [`kache-hash`](https://crates.io/crates/kache-hash) wraps them in a safe interface, which most programs want instead.

The bindings mirror the header `kache_hash/kache_capi.h`, which documents each function: concurrent sets of 31-mers, and maps from 31-mers to unsigned values of 8, 16, 32 or 64 bits.
The library is C++20, compiled by the [`cc`](https://crates.io/crates/cc) crate from the sources in the package, with no CMake: it takes GCC 11 or Clang 12, or later, as `CXX` or the system's `c++`, and an x86-64 CPU with AVX2.
Its build script exports the header's directory to the build scripts of dependents as `DEP_KACHE_CAPI_INCLUDE`.

## Target CPU

The library compiles for:

- the ISA that `KACHE_HASH_ARCH` names, if set: `native`, `x86-64-v3`, `x86-64-v4`, or `none`, for the compiler's own flags, such as those of `CXXFLAGS`;
- else, the native CPU, if rustc targets it, e.g. with `-C target-cpu=native`;
- else, x86-64-v4 (AVX-512), if rustc's target features include AVX-512's F, BW, CD, DQ and VL, e.g. with `-C target-cpu=x86-64-v4`;
- and x86-64-v3 (AVX2) otherwise.

So a default build runs on any x86-64 CPU with AVX2, and the fastest tables for the CPU that builds them come from

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

The library compiles with optimizations in every Cargo profile, as its tables are unusably slow without them.

## License

BSD-3-Clause, as *k*ache-hash's; its vendored [xxHash](https://github.com/Cyan4973/xxHash) is BSD-2-Clause, as in `vendor/include/kache_hash/xxHash/LICENSE`.
