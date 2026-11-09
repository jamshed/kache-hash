
#ifndef KACHE_HASH_VERSION_HPP
#define KACHE_HASH_VERSION_HPP



// kache-hash's version, which its CMake project reads from here. A release
// bumps it, with the crate's and `CITATION.cff`'s.
#define KACHE_HASH_VERSION_MAJOR 0
#define KACHE_HASH_VERSION_MINOR 7
#define KACHE_HASH_VERSION_PATCH 1

// The version as one number, `10000 * major + 100 * minor + patch`, for
// comparisons in `#if`s.
#define KACHE_HASH_VERSION \
    (KACHE_HASH_VERSION_MAJOR * 10000 + KACHE_HASH_VERSION_MINOR * 100 + KACHE_HASH_VERSION_PATCH)

// The version as a string, such as "0.7.1".
#define KACHE_HASH_STRINGIFY_(x) #x
#define KACHE_HASH_STRINGIFY(x) KACHE_HASH_STRINGIFY_(x)
#define KACHE_HASH_VERSION_STRING \
    KACHE_HASH_STRINGIFY(KACHE_HASH_VERSION_MAJOR) "." KACHE_HASH_STRINGIFY(KACHE_HASH_VERSION_MINOR) "." \
    KACHE_HASH_STRINGIFY(KACHE_HASH_VERSION_PATCH)



#endif
