
#ifndef KACHE_HASH_XXHASH_HPP
#define KACHE_HASH_XXHASH_HPP



// xxHash, from its copy in `xxHash/`, with the definitions it needs set here
// rather than by the build. Its functions are inlined into each translation
// unit, as `static` ones prefixed with `XXH_INLINE_`, so that they don't clash
// with another copy of xxHash linked in.


// NB: xxHash reads its inputs with `memcpy` (method 0) here, rather than with
// its default for GCC and Clang (method 1), which breaks strict aliasing: GCC
// 11 at -O3 with AVX2 miscompiles the l-mers' hashes with it. A header that
// compiled xxHash into the translation unit before must have used method 0 too.
#if defined(XXH_IMPLEM_13a8737387) && defined(XXH_FORCE_MEMORY_ACCESS) && XXH_FORCE_MEMORY_ACCESS != 0
#error "kache-hash needs XXH_FORCE_MEMORY_ACCESS 0: include it before xxhash.h, or define that as 0."
#endif

#undef XXH_FORCE_MEMORY_ACCESS
#define XXH_FORCE_MEMORY_ACCESS 0

#ifndef XXH_INLINE_ALL
#define XXH_INLINE_ALL
#endif

#include "xxHash/xxhash.h"



#endif
