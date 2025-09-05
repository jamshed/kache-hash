
#ifndef UTILITY_HPP
#define UTILITY_HPP



#include <cstdint>
#include <cstddef>
#include <cstdlib>
#include <cassert>


// Returns pointer to a memory-allocation for `size` elements of type `T_`.
template <typename T_>
inline T_* allocate(const std::size_t size)
{
    return static_cast<T_*>(std::malloc(size * sizeof(T_)));
}


// Returns pointer to a memory-allocation for `size` elements of type `T_`,
// whose alignment is specified is `alignment`.
template <typename T_>
inline T_* aligned_allocate(const std::size_t size, const std::size_t alignment = 8)
{
    return static_cast<T_*>(std::aligned_alloc(alignment, size * sizeof(T_)));
}


// Returns pointer to a memory-reallocation for `size` elements of type `T_`.
template <typename T_>
inline T_* reallocate(T_* const ptr, const std::size_t size)
{
    return static_cast<T_*>(std::realloc(ptr, size * sizeof(T_)));
}


// Deallocates the pointer `ptr`, allocated with `allocate`.
template <typename T_>
inline void deallocate(T_* const ptr)
{
    std::free(ptr);
}


// Returns the smallest power of 2 at least as large as `x`. `x` must be in
// `[1, 2^63]`.
inline constexpr std::size_t ceil_pow_2(std::size_t x)
{
    assert(x > 0 && x <= (1lu << 63));

    if((x & (x - 1)) == 0)
        return x;

    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);
    x |= (x >> 32);

    return x + 1;
}


// Returns the bits in the indices `[l, h]` (0-indexed) of `x`.
inline constexpr uint64_t bits(uint64_t x, uint64_t l, uint64_t h)
{
    assert(l <= h && h < 64);
    return (x << (63 - h)) >> ((63 - h) + l);
}



#endif
