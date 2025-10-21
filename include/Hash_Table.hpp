
#ifndef HASH_TABLE_HPP
#define HASH_TABLE_HPP



#include "utility.hpp"

#include "xxHash/xxhash.h"

#include <cstdint>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <bit>
#include <fstream>


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
    double lf;  // Maximum load factor.

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

    // Returns the capacity of the set.
    auto capacity() const { return cap; }

    // Clears the set.
    void clear();

    // Inserts `key` into the set and returns `true` iff the insertion succeeds,
    // `key` was not present in the set beforehand.
    bool insert(T_key_ key);

    // Returns `true` iff the set contains `key`.
    bool contains(T_key_ key) const;

    // Serializes the hash set to the stream `os`.
    void serialize(std::ofstream& os) const;

    // Deserializes the hash set from the stream `is`.
    void deserialize(std::ifstream& is);
};


template <typename T_key_, typename T_val_, typename T_hash_ = Hash<T_key_>>
class Hash_Map
{
private:

    static constexpr double lf_default = 0.8;   // Default maximum load factor.
    const double lf;    // Maximum load factor.

    std::size_t cap;    // Capacity of the map.
    std::size_t resize_th;  // Threshold capacity of when to resize the table.

    struct kv_t
    {
        T_key_ key;
        T_val_ val;

        kv_t(const T_key_ key, const T_val_ val): key(key), val(val)
        {}
    }* KV;  // Flat collection of key-value pairs.

    std::size_t sz; // Number of elements in the map.

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

    // Constructs a hash map with an initial capacity `cap` and maximum load
    // factor `lf`.
    Hash_Map(std::size_t cap = 4096, double lf = lf_default);

    ~Hash_Map() { deallocate(KV); }

    Hash_Map(const Hash_Map&) = delete;
    Hash_Map(Hash_Map&&) = delete;
    Hash_Map& operator=(const Hash_Map&) = delete;
    Hash_Map& operator=(Hash_Map&&) = delete;

    // Returns the size of the map.
    auto size() const { return sz; }

    // Returns the capacity of the map.
    auto capacity() const { return cap; }

    // Clears the map.
    void clear();

    // Inserts the key `key` into the map with the associated value `val` only
    // if `key` was not present in the map beforehand. Returns `true` iff the
    // insertion succeeds.
    bool insert(T_key_ key, T_val_ val);

    // Inserts the key `key` into the map with the associated value `T_val()`
    // only if `key` was not present in the map beforehand. Returns a reference
    // to the associated value of `key` in the map.
    T_val_& operator[](T_key_ key);

    // Inserts the key `key`  into the map with the associated value `val` only
    // if `key` was not present in the map beforehand. Otherwise transforms the
    // associated value by applying the function `F` on it. Returns the value
    // of `key` after this operation.
    template <typename T_f_> T_val_ upsert(T_key_ key, T_val_ val, T_f_&& F);

    class Iterator;

    // Returns an iterator to the beginning of the map.
    Iterator begin() const { return Iterator(*this, 0); };

    // Returns an iterator to the end (exclusive) of the map.
    Iterator end() const { return Iterator(*this, cap); };
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
inline void Hash_Set<T_key_, T_hash_>::serialize(std::ofstream& os) const
{
    ::serialize(lf, os);
    ::serialize(cap, os);
    ::serialize(resize_th, os);

    ::serialize(K, cap, os);
    ::serialize(sz, os);
}


template <typename T_key_, typename T_hash_>
inline void Hash_Set<T_key_, T_hash_>::deserialize(std::ifstream& is)
{
    ::deserialize(lf, is);
    ::deserialize(cap, is);
    ::deserialize(resize_th, is);

    if(K)
        deallocate(K);

    K = allocate<T_key_>(cap);
    ::deserialize(K, cap, is);
    ::deserialize(sz, is);
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


template <typename T_key_, typename T_val_, typename T_hash_>
Hash_Map<T_key_, T_val_, T_hash_>::Hash_Map(const std::size_t cap, const double lf):
      lf(lf)
    , cap(ceil_pow_2(cap))
    , resize_th(cap * lf)
    , KV(allocate<kv_t>(cap))
    , hash()
{
    clear();
}


template <typename T_key_, typename T_val_, typename T_hash_>
inline void Hash_Map<T_key_, T_val_, T_hash_>::clear()
{
    std::memset(static_cast<void*>(KV), -1, cap * sizeof(kv_t));
    sz = 0;
}


template <typename T_key_, typename T_val_, typename T_hash_>
inline bool Hash_Map<T_key_, T_val_, T_hash_>::insert(const T_key_ key, const T_val_ val)
{
    for(std::size_t i = map_to_idx(hash(key)); ; i = map_to_idx(i + 1))
        if(KV[i].key == empty)
        {
            KV[i] = kv_t(key, val);
            if(++sz >= resize_th)
                resize();

            return true;
        }
        else if(KV[i].key == key)
            break;

    return false;
}


template <typename T_key_, typename T_val_, typename T_hash_>
inline T_val_& Hash_Map<T_key_, T_val_, T_hash_>::operator[](const T_key_ key)
{
    std::size_t i = map_to_idx(hash(key));
    for(; ; i = map_to_idx(i + 1))
        if(KV[i].key == empty)
        {
            KV[i] = kv_t(key, T_val_());
            if(++sz >= resize_th)
            {
                resize();
                return operator[](key);
            }

            break;
        }
        else if(KV[i].key == key)
            break;

    return KV[i].val;
}


template <typename T_key_, typename T_val_, typename T_hash_>
template <typename T_f_>
inline T_val_ Hash_Map<T_key_, T_val_, T_hash_>::upsert(const T_key_ key, const T_val_ val, T_f_&& F)
{
    for(std::size_t i = map_to_idx(hash(key)); ; i = map_to_idx(i + 1))
        if(KV[i].key == empty)
        {
            KV[i] = KV(key, val);
            if(++sz >= resize_th)
                resize();

            break;
        }
        else if(KV[i].key == key)
            return KV[i].val = F(val);

    return val;
}


template <typename T_key_, typename T_val_, typename T_hash_>
inline void Hash_Map<T_key_, T_val_, T_hash_>::resize()
{
    const auto K_old = KV;
    KV = allocate<kv_t>(cap * 2);
    cap *= 2, resize_th = cap * lf;
    clear();

    for(std::size_t i = 0; i < cap / 2; ++i)
        if(K_old[i].key != empty)
            insert(K_old[i].key, K_old[i].val);

    deallocate(K_old);
}


// Iterator over `Hash_Map`.
template <typename T_key_, typename T_val_, typename T_hash_>
class Hash_Map<T_key_, T_val_, T_hash_>::Iterator
{
    friend class Hash_Map<T_key_, T_val_, T_hash_>;

    typedef Hash_Map<T_key_, T_val_, T_hash_> ht_t;

private:

    const ht_t& ht; // Associated hash map.
    std::size_t off;    // Position of the iterator in the map.

    // Constructs an iterator for the hash map `ht` at offset `off`.
    Iterator(const ht_t& ht, std::size_t off): ht(ht), off(off)
    {
        assert(off == 0 || off == ht.cap);
        if(off == 0 && ht.KV[0].key == ht.empty)
            operator++();
    }

public:

    // Repositions the iterator to the next successive element in the map, or to
    // its end if there are no more elements.
    void operator++()
    {
        if(off == ht.cap)
            return;

        do
        {
            off++;
            if(ht.KV[off].key != ht.empty)
                break;
        }
        while(off < ht.cap);
    }

    const ht_t::kv_t& operator*() { assert(off < ht.cap); return ht.KV[off]; }

    bool operator==(const Iterator& rhs) const { return &ht == &rhs.ht && off == rhs.off; }
};



#endif
