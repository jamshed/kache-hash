
#ifndef KACHE_HASH_ARCH_HPP
#define KACHE_HASH_ARCH_HPP



// The SIMD intrinsics of the target architecture: AVX2 on x86-64, and NEON
// (Advanced SIMD), which every AArch64 core has, on AArch64.
#if defined(__x86_64__) || defined(_M_X64)
#define KACHE_HASH_X86_64 1
#include <immintrin.h>
#elif defined(__aarch64__) || defined(_M_ARM64)
#define KACHE_HASH_AARCH64 1
#include <arm_neon.h>
#else
#error "kache-hash supports x86-64 (with AVX2) and AArch64 (with NEON) only."
#endif


namespace kache_hash
{

// Hints the core that the thread is spinning on a wait, e.g. for a lock.
inline void cpu_relax()
{
#if defined(KACHE_HASH_X86_64)
    _mm_pause();
#else
    // NB: `yield` is a no-op on most AArch64 cores; `isb` delays the spinner for a few cycles, somewhat as `pause`.
    __asm__ __volatile__("isb");
#endif
}

}



#endif
