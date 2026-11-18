# Installation

*k*ache-hash takes a C++20 compiler, GCC 11 or later or Clang 12 or later, and an x86-64 CPU with AVX2.
Its table is header-only.
Its C interface and its counter are compiled libraries, which CMake builds, and the Rust crates build the C interface's library themselves.

## CMake

CMake 3.14 or later builds the libraries and installs them with the headers, a CMake package, and pkg-config files:

```bash
cmake -S . -B build
cmake --build build -j
cmake --install build --prefix <prefix>
```

A top-level build is a Release build unless `CMAKE_BUILD_TYPE` says otherwise, and it builds the tests, which need zlib, and the examples too.

### Options

| Option | Default | Effect |
|---|---|---|
| `KACHE_HASH_ARCH` | `native` | The ISA of the compiled libraries: `native`, the building CPU's; `x86-64-v3` (AVX2) or `x86-64-v4` (AVX-512), for binaries that run on other CPUs, with `x86-64-v3` the default under conda-build and cibuildwheel; or `none`, for the toolchain's flags alone. |
| `BUILD_SHARED_LIBS` | `OFF` | Builds the libraries shared, exporting the symbols of their interfaces only, rather than static. |
| `KACHE_HASH_BUILD_CAPI` | `ON` | Builds the C interface's library, `kache_capi`. |
| `KACHE_HASH_BUILD_COUNTER` | `ON` | Builds the counter's library, `kache_counter`. |
| `KACHE_HASH_INSTALL` | in top-level builds | Generates the install rules. |
| `KACHE_HASH_BUILD_TESTS` | in top-level builds | Builds the tests, for CTest. |
| `KACHE_HASH_BUILD_EXAMPLES` | in top-level builds | Builds the examples, which CTest runs with the tests. |
| `KACHE_HASH_BUILD_BENCHMARKS` | `OFF` | Builds the benchmarks and the tools, fetching the baselines; see [Benchmarks](benchmarks.md). |
| `KACHE_HASH_TEST_THREADS` | all CPUs | The threads of the concurrent tests. |
| `DEBUG_INFO`, `ASAN`, `UBSAN`, `TSAN` | `OFF` | Development builds: debug information, and the sanitizers of addresses, undefined behavior (with Clang), and threads, which goes alone. |

Branches of the compiled code are kept from crossing or ending on 32-byte boundaries, for the JCC erratum of Skylake-derived cores, such as Cascade Lake: that gains 1–4% in builds there, and up to 9% in negative queries.

### Using an installed *k*ache-hash

A project finds the installed package, with the components it needs, and links its targets:

```cmake
find_package(kache_hash 0.8 REQUIRED COMPONENTS capi counter)
target_link_libraries(app PRIVATE kache_hash::kache_hash)   # The header-only table.
target_link_libraries(tool PRIVATE kache_hash::kache_capi)  # The C interface.
```

`kache_hash::kache_hash` carries the include path, C++20, and threads, while the optimization and ISA flags, AVX2 at least, stay the project's own.
The libraries' interfaces are plain C and C++, without the table's headers, so a project linking them keeps its own flags.
The package's version matches releases of the same minor version, as *k*ache-hash is at 0.x.

`FetchContent` takes an installed *k*ache-hash before fetching its sources:

```cmake
include(FetchContent)
FetchContent_Declare(kache_hash GIT_REPOSITORY https://github.com/jamshed/kache-hash.git GIT_TAG v0.8.0 FIND_PACKAGE_ARGS)
FetchContent_MakeAvailable(kache_hash)
```

A project that builds *k*ache-hash's sources this way gets neither its tests nor its examples, and its install rules only with `KACHE_HASH_INSTALL`.

### pkg-config

`kache_hash.pc` gives the flags of the headers, and `kache_capi.pc` those of the C interface, with the C++ standard library that a static `kache_capi` needs:

```bash
cc -std=c11 app.c $(pkg-config --cflags --libs kache_capi)
c++ -std=c++20 -mavx2 -O3 app.cpp $(pkg-config --cflags --libs kache_hash)
```

## Conda

```bash
conda install -c conda-forge -c bioconda kache-hash
```

Bioconda's package has the headers, the static and shared libraries for x86-64-v3, the CMake package, and the pkg-config files.

## Cargo

```bash
cargo add kache-hash
```

The crate builds the C interface's library through the crate `kache-hash-sys`, with the `cc` crate rather than CMake; [the Rust guide](rust.md) tells which ISA it targets.

## The headers alone

Copy `include/kache_hash` into a directory on the project's include path, and include `<kache_hash/kache_hash.hpp>`, or `<kache_hash/Streaming_Kmer_Hash_Table.hpp>` for the table alone.
The headers need only C++20 and AVX2, e.g. `-std=c++20 -mavx2`, and no definitions.
They assume 64-byte cache lines, which `KACHE_HASH_CACHE_LINE_SIZE` changes.
