# Supported platforms

- **CPUs.** x86-64 with AVX2, i.e. Intel's from Haswell and AMD's from Excavator and Zen. A build for AVX-512, i.e. for `x86-64-v4`, or `native` on a CPU with it, uses it too, in its SIMD lanes and checksums, and runs only on CPUs with it.
- **Systems.** Linux, with GCC 11 or later or Clang 12 or later; and macOS on Intel, with Apple Clang. CI tests both, and AVX-512 under Intel's SDE.
- **Apple silicon.** Under Rosetta 2, which runs AVX2 code from macOS 15: an x86-64 build, e.g. CMake's with `-DCMAKE_OSX_ARCHITECTURES=x86_64`, whose `native` targets Haswell there, or `cargo build --target x86_64-apple-darwin`, runs there. aarch64 itself, with NEON, is planned for v0.9.0.

The headers include `<immintrin.h>` and stop with an error without AVX2; the compiled libraries take their ISA from `KACHE_HASH_ARCH`, as [Installation](install.md) describes.
