
#ifndef KACHE_HASH_CONCURRENT_HASH_TABLE_HPP
#define KACHE_HASH_CONCURRENT_HASH_TABLE_HPP



#include "Spin_Lock.hpp"
#include "utility.hpp"
#include "xxHash/xxhash.h"

#include <cstddef>
#include <cstring>
#include <optional>
#include <cstdlib>
#include <algorithm>
#include <cmath>


namespace kache_hash
{


// Hasher for objects of type `T_`.
template <typename T_>
struct Hash
{
    uint64_t operator()(const T_& key) const { return XXH3_64bits(&key, sizeof(key)); }
};


// =============================================================================
template <typename T_key_, typename T_val_, typename T_hasher_>
class Concurrent_Hash_Table
{
    static_assert(sizeof(T_key_) <= 8, "Key for concurrent hash table must have size at most 8 bytes for atomic-read guarantees to hold.");

    class Iterator;
    friend class Iterator;

private:

    // A slot. Its key is stored complemented, so that zeroed memory, fresh from
    // the operating system, holds empty slots.
    struct Key_Val_Pair
    {
        T_key_ key;
        [[no_unique_address]] T_val_ val;   // An empty value, such as a set's, takes no room.
    };

    // Returns the bitwise complement of the key `key`: a key's stored form.
    static T_key_ stored(T_key_ key);

    static constexpr double lf_default_ = 0.8;  // Default maximum load-factor supported.

    T_key_ empty_key_;  // The stored form of empty slots: all 0-bits, the stored form of the all-1-bits key.

    const T_hasher_ hash;   // The hasher object.

    const std::size_t capacity_;    // True capacity of the table; adjusted to be a power of 2.
    const std::size_t idx_wrapper_mask; // Bitmask to wrap indexing into the table.

    Key_Val_Pair* const T;  // The flat table of key-value collection.

    Spin_Lock* const lock;  // Locks for atomic storing of key-value entries; zeroed memory holds them unlocked.
    // TODO: experiment moving the spin-lock to the key-value structure for seemingly better cache-behavior.

    // Maps the hash value `h` to an index into the table.
    std::size_t hash_to_idx(std::size_t h) const { return h & idx_wrapper_mask; }

    // Returns the next (wrapped) index for `i`.
    std::size_t next_index(std::size_t i) const { return hash_to_idx(i + 1); }


public:

    // Constructs a concurrent hash table to support upto `max_n` elements, with
    // a maximum load-factor of `load_factor`. The object `hasher` is used to
    // hash the keys in the table. Mapping the table takes `thread_c` threads.
    Concurrent_Hash_Table(std::size_t max_n, std::size_t thread_c = 1, double load_factor = lf_default_, T_hasher_ hasher = T_hasher_());

    ~Concurrent_Hash_Table();

    // Returns the capacity of the hash table.
    std::size_t capacity() const { return capacity_; }

    // Returns the default maximum load-factor.
    static constexpr double lf_default() { return lf_default_; }

    // Returns the size of the table. It costs linear in the size of the table.
    std::size_t size() const;

    // Inserts the key `key` with value `val` into the table. Returns `nullptr`
    // iff the insertion succeeds, i.e. the key was absent in the table. Returns
    // the address of the associated value otherwise. Sets `full` iff the key
    // is absent and the table has no slot left for it; then nothing changes
    // and `nullptr` is returned.
    const T_val_* insert(T_key_ key, T_val_ val, bool& full);

    // Inserts the key `key` with value `val` into the table if it is absent and
    // returns `null_opt`. Otherwise updates the associated value `v` to `f(v)`
    // and returns the old value. Sets `full` iff the key is absent and the
    // table has no slot left for it; then nothing changes.
    template <typename F_> std::optional<T_val_> upsert(T_key_ key, T_val_ val, const F_& f, bool& full);

    // Searches for `key` in the table and returns the address of the value
    // associated to it iff it is found. Returns `nullptr` otherwise.
    const T_val_* find(T_key_ key) const;

    // Prefetches the first slot probed for `key`, and its lock iff `write_`,
    // with write intent iff `write_`.
    template <bool write_> void prefetch(T_key_ key) const;

    // Returns an iterator for the key-value pairs in the table that belongs to
    // a group of `it_count` iterators and has an ID `it_id` in the group.
    Iterator iterator(const std::size_t it_count, const std::size_t it_id) { return Iterator(*this, it_count, it_id); }

    // Returns the memory of the range of slots of the iterator with ID `it_id`
    // in a group of `it_count` to the operating system, where supported, which
    // empties them. For a table about to be discarded.
    void release(std::size_t it_count, std::size_t it_id);
};


// =============================================================================
template <typename T_key_, typename T_val_, typename T_hasher_>
class Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::Iterator
{
    friend class Concurrent_Hash_Table;

private:

    Concurrent_Hash_Table& M;   // The hash table to iterate on.
    std::size_t idx;    // Current slot-index the iterator is in.
    std::size_t end;    // End-index of the slot-range for the iterator.


    // Constructs an iterator for the hash table `M`, that would belong to a
    // group of `it_count` iterators and has an ID `it_id` in the group.
    Iterator(Concurrent_Hash_Table& M, std::size_t it_count, std::size_t it_id);


public:

    // Moves the iterator to the next key in the table. Iff some key is found
    // within its remaining range, the key and the associated value are put at
    // `key` and `val` respectively, and returns true.
    bool next(T_key_& key, T_val_& val);
};


template <typename T_key_, typename T_val_, typename T_hasher_>
inline Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::Concurrent_Hash_Table(const std::size_t max_n, const std::size_t thread_c, const double load_factor, const T_hasher_ hasher):
      empty_key_()
    , hash(hasher)
    , capacity_(ceil_pow_2(static_cast<std::size_t>(std::ceil(max_n / load_factor))))
    , idx_wrapper_mask(capacity_ - 1)
    , T(allocate_huge<Key_Val_Pair>(capacity_))
    , lock(allocate_huge<Spin_Lock>(capacity_))
{
    std::memset(reinterpret_cast<char*>(&empty_key_), 0, sizeof(empty_key_));   // The fresh slots are empty already.

    // Map the table up front, in parallel; else probes would first read its pages, mapping the shared huge zero page,
    // which the first writes then split into small pages.
    parallel_prefault(T, capacity_ * sizeof(Key_Val_Pair), thread_c);
    parallel_prefault(lock, capacity_ * sizeof(Spin_Lock), thread_c);
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::~Concurrent_Hash_Table()
{
    deallocate_huge(T, capacity_), deallocate_huge(lock, capacity_);
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline T_key_ Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::stored(const T_key_ key)
{
    T_key_ x;
    auto* const b = reinterpret_cast<unsigned char*>(&x);
    std::memcpy(b, &key, sizeof(key));
    for(std::size_t i = 0; i < sizeof(key); ++i)
        b[i] = ~b[i];
    return x;
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline std::size_t Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::size() const
{
    std::size_t sz = 0;
    for(std::size_t i = 0; i < capacity_; ++i)
        sz += (T[i].key != empty_key_);

    return sz;
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline const T_val_* Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::insert(const T_key_ key, const T_val_ val, bool& full)
{
    bool success = false;
    full = false;
    const auto key_st = stored(key);

    std::size_t i = hash_to_idx(hash(key));
    for(std::size_t tried_slots = 0; tried_slots < capacity_; ++tried_slots, i = next_index(i))
    {
        if(T[i].key == empty_key_)
        {
            lock[i].lock();
            if(T[i].key == empty_key_)
                T[i].key = key_st, T[i].val = val,
                success = true;
            lock[i].unlock();

            if(success)
                return nullptr;
        }

        if(T[i].key == key_st)
            return &T[i].val;
    }

    full = true;
    return nullptr;
}


template <typename T_key_, typename T_val_, typename T_hasher_>
template <typename F_>
inline std::optional<T_val_> Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::upsert(const T_key_ key, const T_val_ val, const F_& f, bool& full)
{
    bool success = false;
    full = false;
    const auto key_st = stored(key);

    std::size_t i = hash_to_idx(hash(key));
    for(std::size_t tried_slots = 0; tried_slots < capacity_; ++tried_slots, i = next_index(i))
    {
        if(T[i].key == empty_key_)
        {
            lock[i].lock();
            if(T[i].key == empty_key_)
                T[i].key = key_st, T[i].val = val,
                success = true;
            lock[i].unlock();

            if(success)
                return std::nullopt;
        }

        if(T[i].key == key_st)
        {
            lock[i].lock();
            const auto old_val = T[i].val;
            T[i].val = f(old_val);
            lock[i].unlock();

            return old_val;
        }
    }

    full = true;
    return std::nullopt;
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline const T_val_* Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::find(const T_key_ key) const
{
    const auto key_st = stored(key);
    const T_val_* val_ptr = nullptr;
    std::size_t i = hash_to_idx(hash(key));
    for(std::size_t tried_slots = 0; tried_slots < capacity_; ++tried_slots, i = next_index(i))
    {
        if(T[i].key == key_st)
        {
            val_ptr = &T[i].val;    // NB: The caller reads the value after the return anyway, so locking the slot here guards nothing.
            break;
        }
        else if(T[i].key == empty_key_)
            break;
    }

    return val_ptr;
}


template <typename T_key_, typename T_val_, typename T_hasher_>
template <bool write_>
inline void Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::prefetch(const T_key_ key) const
{
    const auto i = hash_to_idx(hash(key));
    __builtin_prefetch(&T[i], write_, 3);
    if constexpr(write_)
        __builtin_prefetch(&lock[i], write_, 3);
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline void Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::release(const std::size_t it_count, const std::size_t it_id)
{
    const auto range_sz = (capacity_ + it_count - 1) / it_count;
    const auto b = std::min(it_id * range_sz, capacity_), e = std::min((it_id + 1) * range_sz, capacity_);
    release_huge(T + b, e - b), release_huge(lock + b, e - b);
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::Iterator::Iterator(Concurrent_Hash_Table& M, std::size_t it_count, std::size_t it_id):
      M(M)
{
    const auto range_sz = (M.capacity_ + it_count - 1) / it_count;
    idx = it_id * range_sz;
    end = std::min((it_id + 1) * range_sz, M.capacity_);
}


template <typename T_key_, typename T_val_, typename T_hasher_>
inline bool Concurrent_Hash_Table<T_key_, T_val_, T_hasher_>::Iterator::next(T_key_& key, T_val_& val)
{
    for(; idx < end; idx++)
        if(M.T[idx].key != M.empty_key_)
        {
            key = stored(M.T[idx].key), val = M.T[idx].val;
            idx++;
            return true;
        }

    return false;
}


}



#endif