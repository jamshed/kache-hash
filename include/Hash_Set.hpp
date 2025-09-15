
#ifndef HASH_SET_HPP
#define HASH_SET_HPP



#include "utility.hpp"

#include "xxHash/xxhash.h"

#include <cstdint>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <bit>


// Hasher for objects of type `T_`.
template <typename T_>
struct Hash
{
    uint64_t operator()(const T_& key) const { return XXH3_64bits(&key, sizeof(key)); }
};


template <typename T_key_, typename T_hash_ = Hash<T_key_>>
class Hash_Set
{
private:

    static constexpr double lf_default = 0.8;   // Default maximum load factor.
    const double lf;    // Maximum load factor.

    std::size_t cap;    // Capacity of the set.
    std::size_t resize_th;  // Threshold capacity of when to resize the table.

    T_key_* K;  // Flat collection of keys.
    std::size_t sz; // Number of keys in the set.

    T_hash_ hash;   // Hasher.

    static constexpr T_key_ empty = // The empty key.
        []()
        {
            constexpr std::size_t sz = sizeof(T_key_);
            if constexpr(sz == 1)       return std::bit_cast<T_key_>(int8_t(-1));
            else if constexpr(sz == 2)  return std::bit_cast<T_key_>(int16_t(-1));
            else if constexpr(sz == 4)  return std::bit_cast<T_key_>(int32_t(-1));
            else if constexpr(sz == 8)  return std::bit_cast<T_key_>(int64_t(-1));
            else if constexpr(sz == 16) return std::bit_cast<T_key_>(__int128_t(-1));
            else std::exit(EXIT_FAILURE);
        }();


    // Maps (i.e. wraps) the value `h` to the appropriate index into the table.
    std::size_t map_to_idx(const uint64_t h) const { return h & (cap - 1); }

    // Expands the table geometrically.
    void resize();


public:

    // Constructs a hash set with an initial capacity `cap` and maximum load
    // factor `lf`.
    Hash_Set(std::size_t cap = 4096, double lf = lf_default);

    ~Hash_Set() { deallocate(K); }

    Hash_Set(const Hash_Set&) = delete;
    Hash_Set(Hash_Set&&) = delete;
    Hash_Set& operator=(const Hash_Set&) = delete;
    Hash_Set& operator=(Hash_Set&&) = delete;

    // Returns the size of the set.
    auto size() const { return sz; }

    // Clears the set.
    void clear();

    // Inserts `key` into the set and returns `true` iff the insertion succeeds,
    // `key` was not present in the set beforehand.
    bool insert(T_key_ key);

    // Returns `true` iff the set contains `key`.
    bool contains(T_key_ key) const;
};


template <typename T_key_, typename T_hash_>
Hash_Set<T_key_, T_hash_>::Hash_Set(const std::size_t cap, const double lf):
      lf(lf)
    , cap(ceil_pow_2(cap))
    , resize_th(cap * lf)
    , K(allocate<T_key_>(cap))
    , hash()
{
    clear();
}


template <typename T_key_, typename T_hash_>
inline void Hash_Set<T_key_, T_hash_>::clear()
{
    std::memset(static_cast<void*>(K), -1, cap * sizeof(T_key_));
    sz = 0;
}


template <typename T_key_, typename T_hash_>
inline bool Hash_Set<T_key_, T_hash_>::insert(const T_key_ key)
{
    for(std::size_t i = map_to_idx(hash(key)); ; i = map_to_idx(i + 1))
        if(K[i] == empty)
        {
            K[i] = key;
            if(++sz >= resize_th)
                resize();

            return true;
        }
        else if(K[i] == key)
            break;

    return false;
}


template <typename T_key_, typename T_hash_>
inline bool Hash_Set<T_key_, T_hash_>::contains(const T_key_ key) const
{
    for(std::size_t i = map_to_idx(hash(key)); ; i = map_to_idx(i + 1))
        if(K[i] == empty)
            break;
        else if(K[i] == key)
            return true;

    return false;
}


template <typename T_key_, typename T_hash_>
inline void Hash_Set<T_key_, T_hash_>::resize()
{
    const auto K_old = K;
    K = allocate<T_key_>(cap * 2);
    cap *= 2, resize_th = cap * lf;
    clear();

    for(std::size_t i = 0; i < cap / 2; ++i)
        if(K_old[i] != empty)
            insert(K_old[i]);

    deallocate(K_old);
}



#endif
