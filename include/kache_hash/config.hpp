
#ifndef KACHE_HASH_CONFIG_HPP
#define KACHE_HASH_CONFIG_HPP



// What the headers take from the compiler, with their defaults: they build with
// C++20 and AVX2 alone, e.g. with `-std=c++20 -mavx2`, and no definitions.


#ifndef __AVX2__
#error "kache-hash needs AVX2: compile with -mavx2, -march=x86-64-v3, or -march=native on a CPU with AVX2."
#endif


// The bytes of a cache line, which the tables align and pad their parts to: 64
// by default, as on x86-64 CPUs. `KACHE_HASH_CACHE_LINE_SIZE` sets it, as does
// `L1_CACHE_LINE_SIZE`, which builds defined before v0.8.0.
#if !defined(KACHE_HASH_CACHE_LINE_SIZE) && defined(L1_CACHE_LINE_SIZE)
#define KACHE_HASH_CACHE_LINE_SIZE L1_CACHE_LINE_SIZE
#elif !defined(KACHE_HASH_CACHE_LINE_SIZE)
#define KACHE_HASH_CACHE_LINE_SIZE 64
#endif



#endif
