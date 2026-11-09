
#ifndef KACHE_HASH_UTILITY_HPP
#define KACHE_HASH_UTILITY_HPP



#include "config.hpp"

#include <cstdint>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <thread>
#include <vector>
#include <algorithm>
#include <cassert>

#ifdef __linux__
#include <sys/mman.h>
#endif


namespace kache_hash
{


constexpr std::size_t huge_page_sz = std::size_t(1) << 21;  // Size of the transparent huge pages targeted, 2 MB.


// Returns the number of bytes `allocate_huge` maps for `size` elements of type
// `T_`.
template <typename T_>
inline std::size_t huge_bytes(const std::size_t size)
{
    return ((size * sizeof(T_) + huge_page_sz - 1) / huge_page_sz) * huge_page_sz;
}


// Returns pointer to a memory-allocation for `size` elements of type `T_`,
// aligned to huge pages and advised to be backed by them where supported. The
// memory is zeroed. Release it with `deallocate_huge`.
template <typename T_>
inline T_* allocate_huge(const std::size_t size)
{
    const auto bytes = huge_bytes<T_>(size);
#ifdef __linux__
    // Over-map by a huge page, then trim to an aligned range.
    auto* const raw = static_cast<char*>(mmap(nullptr, bytes + huge_page_sz, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0));
    if(raw == MAP_FAILED)
        return nullptr;

    auto* const p = reinterpret_cast<char*>((reinterpret_cast<uintptr_t>(raw) + huge_page_sz - 1) & ~(huge_page_sz - 1));
    if(p > raw)
        munmap(raw, p - raw);
    munmap(p + bytes, (raw + bytes + huge_page_sz) - (p + bytes));
    madvise(p, bytes, MADV_HUGEPAGE);
    return reinterpret_cast<T_*>(p);
#else
    auto* const p = static_cast<T_*>(std::aligned_alloc(huge_page_sz, bytes));
    if(p)
        std::memset(static_cast<void*>(p), 0, bytes);
    return p;
#endif
}


// Deallocates the pointer `ptr`, allocated with `allocate_huge` for `size`
// elements of type `T_`.
template <typename T_>
inline void deallocate_huge(T_* const ptr, const std::size_t size)
{
    if(!ptr)
        return;
#ifdef __linux__
    munmap(static_cast<void*>(ptr), huge_bytes<T_>(size));
#else
    (void)size;
    std::free(ptr);
#endif
}


// Grows the allocation of `allocate_huge` at `ptr` from `old_size` to
// `new_size` elements of type `T_`, keeping its contents, and returns its
// possibly moved address; the elements added are zeroed. On Linux the pages
// move to an aligned range of their own without copying, huge pages
// included.
template <typename T_>
inline T_* grow_huge(T_* const ptr, const std::size_t old_size, const std::size_t new_size)
{
    const auto old_bytes = huge_bytes<T_>(old_size), new_bytes = huge_bytes<T_>(new_size);
#ifdef __linux__
    // Reserve an aligned range, then move the pages there, extended.
    auto* const raw = static_cast<char*>(mmap(nullptr, new_bytes + huge_page_sz, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0));
    if(raw == MAP_FAILED)
        return nullptr;

    auto* const dst = reinterpret_cast<char*>((reinterpret_cast<uintptr_t>(raw) + huge_page_sz - 1) & ~(huge_page_sz - 1));
    if(dst > raw)
        munmap(raw, dst - raw);
    munmap(dst + new_bytes, (raw + new_bytes + huge_page_sz) - (dst + new_bytes));

    void* const p = mremap(static_cast<void*>(ptr), old_bytes, new_bytes, MREMAP_MAYMOVE | MREMAP_FIXED, dst);
    if(p == MAP_FAILED)
        return nullptr;

    madvise(p, new_bytes, MADV_HUGEPAGE);
    return static_cast<T_*>(p);
#else
    (void)new_bytes;
    auto* const p = allocate_huge<T_>(new_size);
    if(p)
        std::memcpy(static_cast<void*>(p), static_cast<const void*>(ptr), old_bytes);
    deallocate_huge(ptr, old_size);
    return p;
#endif
}


// Returns the memory of the `size` elements of type `T_` at `ptr`, within an
// allocation of `allocate_huge`, to the operating system, where supported:
// they read as zeros thereafter. The whole huge pages within the range are
// released.
template <typename T_>
inline void release_huge(T_* const ptr, const std::size_t size)
{
#ifdef __linux__
    const auto b = (reinterpret_cast<uintptr_t>(ptr) + huge_page_sz - 1) & ~(huge_page_sz - 1);
    const auto e = (reinterpret_cast<uintptr_t>(ptr + size)) & ~(huge_page_sz - 1);
    if(b < e)
        madvise(reinterpret_cast<void*>(b), e - b, MADV_DONTNEED);
#else
    (void)ptr, (void)size;
#endif
}


// Maps each page of the `bytes` bytes of fresh, zeroed memory at `ptr` up
// front, in parallel with `thread_c` threads, by writing a zero to it. A write
// is the first touch on purpose: a read first would map the shared huge zero
// page, which a later write splits into small pages.
inline void parallel_prefault(void* const ptr, const std::size_t bytes, const std::size_t thread_c)
{
    constexpr std::size_t page_sz = 4096;
    constexpr std::size_t min_share = std::size_t(1) << 24;  // Smallest share worth a thread, 16 MB.
    const auto t_c = std::max<std::size_t>(1, std::min(thread_c, bytes / min_share));
    const auto touch = [=](const std::size_t s, const std::size_t e)
    {
        auto* const p = static_cast<volatile char*>(ptr);
        for(std::size_t off = s; off < e; off += page_sz)
            p[off] = 0;
    };

    if(t_c == 1)
        return touch(0, bytes);

    std::vector<std::thread> W;
    W.reserve(t_c);
    for(std::size_t t = 0; t < t_c; ++t)
        W.emplace_back(touch, (bytes / huge_page_sz * t / t_c) * huge_page_sz, (t + 1 == t_c ? bytes : (bytes / huge_page_sz * (t + 1) / t_c) * huge_page_sz));

    std::for_each(W.begin(), W.end(), [](auto& w){ w.join(); });
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


// Fields for profiling time.
constexpr static auto now = std::chrono::high_resolution_clock::now;
constexpr static auto duration = [](const std::chrono::nanoseconds& d) { return std::chrono::duration_cast<std::chrono::duration<double>>(d).count(); };


// Wrapper class for a data element of type `T_` to ensure that in a linear
// collection of `T_`'s, each element is aligned to a cache-line boundary.
template <typename T_>
class alignas(KACHE_HASH_CACHE_LINE_SIZE)
    Padded
{
private:

    T_ data_;


public:

    Padded()
    {}

    Padded(const T_& data):
      data_(data)
    {}

    Padded(T_&& data):
        data_(std::move(data))
    {}

    T_& unwrap() { return data_; }

    const T_& unwrap() const { return data_; }
};


}



#endif
