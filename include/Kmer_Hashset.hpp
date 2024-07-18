
#ifndef KMER_HASHSET_HPP
#define KMER_HASHSET_HPP



#include "Kmer.hpp"
#include "utility.hpp"

#include <cstdint>
#include <cstddef>
#include <limits>
#include <cstring>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cassert>


// A fixed-sized hashset for `k`-mers.
template <uint16_t k>
class Kmer_Hashset
{
private:

    static constexpr double lf_default = 0.75;  // Default maximum load-factor.

    const std::size_t capacity_;    // True capacity of the set; adjusted to be a power of 2.
    const std::size_t wrapper_mask; // Bitmask to wrap indexing into the table of the set.

    Kmer<k>* const K;   // Key collection.
    uint8_t* const T;   // Timestamp of the table version where each entry belongs to.

    std::size_t sz; // Size of the hashset.

    uint8_t cur_stamp;  // Current timestamp, or the version of the table. It need not be just 1 byte.

    std::vector<Kmer<k>> U; // Buffer of update-entries.
    static constexpr std::size_t buf_cap = 1 * 1024 * 1024 / (sizeof(Kmer<k>)); // 1MB worth of update entries.


    // Maps the hash value `h` to an index into the table.
    auto hash_to_idx(const std::size_t h) const { return h & wrapper_mask; }

    // Returns the next (wrapped) index to `i`.
    auto next_idx(const std::size_t idx) const { return hash_to_idx(idx + 1); }

    // Attemps to insert the `k`-mer `key` onward from the index `idx`. Returns
    // `true` iff the insertion was successful, i.e. when the `k`-mer was not
    // present in the set.
    bool insert_at(const Kmer<k>& key, std::size_t idx);


public:

    // Constructs a hashset to support upto `max_n` `k`-mers, with a maximum
    // load-factor of `lf`.
    Kmer_Hashset(std::size_t max_n, double lf = lf_default);

    Kmer_Hashset(const Kmer_Hashset&) = delete;
    Kmer_Hashset(Kmer_Hashset&&);
    Kmer_Hashset& operator=(const Kmer_Hashset&) = delete;
    Kmer_Hashset& operator=(Kmer_Hashset&&) = delete;

    ~Kmer_Hashset() { deallocate(K), deallocate(T); }

    // Returns the size of the hashset.
    auto size() const { return sz; }

    // Returns the true capacity of the hashset.
    auto capacity() const { return capacity_; }

    // Clears the hashset.
    void clear();

    // Inserts the `k`-mer `key` to the hashset. Returns `true` iff the
    // insertion was successful, i.e. when the `k`-mer was not present in the
    // set.
    bool insert(const Kmer<k>& key);

    // Adds the `k`-mer `key` to the hashset buffer for insertion at some time.
    // Returns `true` iff this call triggers a buffer flush.
    bool insert_lazy(const Kmer<k>& key);

    // Completes and flushes all the pending insertions of the hashset.
    void flush_buffer();
};


template <uint16_t k>
inline Kmer_Hashset<k>::Kmer_Hashset(const std::size_t max_n, const double lf):
      capacity_(ceil_pow_2(static_cast<std::size_t>(std::ceil(max_n / lf))))
    , wrapper_mask(capacity_ - 1)
    , K(allocate<Kmer<k>>(capacity_))
    , T(allocate<uint8_t>(capacity_))
    , sz(0)
    , cur_stamp(std::numeric_limits<uint8_t>::max())
{
    assert(capacity_ > 1);

    clear();
}


template <uint16_t k>
inline Kmer_Hashset<k>::Kmer_Hashset(Kmer_Hashset&& rhs):
      capacity_(rhs.capacity_)
    , wrapper_mask(rhs.wrapper_mask)
    , K(rhs.K)
    , T(rhs.T)
    , sz(rhs.sz)
    , cur_stamp(rhs.cur_stamp)
{
    // Moved objects are not really moved in C++.
    const_cast<Kmer<k>*&>(rhs.K) = nullptr; // NOLINT(cppcoreguidelines-pro-type-const-cast)
    const_cast<Kmer<k>*&>(rhs.T) = nullptr; // NOLINT(cppcoreguidelines-pro-type-const-cast)
}


template <uint16_t k>
inline void Kmer_Hashset<k>::clear()
{
    sz = 0;
    cur_stamp = (cur_stamp + 1u) & std::numeric_limits<decltype(cur_stamp)>::max();
    if(cur_stamp == 0)
    {
        std::memset(reinterpret_cast<char*>(T), 0, capacity_ * sizeof(cur_stamp));
        cur_stamp = 1;
    }
}


template <uint16_t k>
inline bool Kmer_Hashset<k>::insert(const Kmer<k>& key)
{
    return insert_at(key, hash_to_idx(key.to_u64()));
}


template <uint16_t k>
inline bool Kmer_Hashset<k>::insert_at(const Kmer<k>& key, const std::size_t idx)
{
    for(std::size_t i = idx; ; i = next_idx(i))
        if(T[i] != cur_stamp)   // Empty slot as it contains an entry from a previous version of the table.
        {
            assert(size() < capacity());

            T[i] = cur_stamp;
            K[i] = key;
            sz++;

            return true;
        }
        else if(K[i] == key)
            return false;

    return false;
}


template <uint16_t k>
inline bool Kmer_Hashset<k>::insert_lazy(const Kmer<k>& key)
{
    U.emplace_back(key);
    if(U.size() == buf_cap)
    {
        flush_buffer();
        return true;
    }

    return false;
}


template <uint16_t k>
inline void Kmer_Hashset<k>::flush_buffer()
{
    constexpr std::size_t batch_size = 64;
    uint64_t H[2][batch_size];

    std::size_t processed = 0;
    std::size_t cur_batch_sz;
    std::size_t next_batch_sz = std::min(batch_size, U.size());
    for(std::size_t i = 0; i < next_batch_sz; ++i)
        H[0][i] = hash_to_idx(U[i].to_u64());

    for(std::size_t b_id = 0; processed < U.size(); ++b_id)
    {
        cur_batch_sz = next_batch_sz;
        const auto next_g_base = processed + cur_batch_sz;
        next_batch_sz = std::min(batch_size, U.size() - next_g_base);
        const auto* const h_cur = &H[b_id & 1][0];
        auto* const h_next = &H[(b_id + 1) & 1][0];

        // Prefetch next batch's initial slots.
        for(std::size_t i = 0; i < next_batch_sz; ++i)
        {
            h_next[i] = hash_to_idx(U[next_g_base + i].to_u64());
            __builtin_prefetch(T + h_next[i], 1, 0);
        }

        // Process current batch.
        for(std::size_t i = 0; i < cur_batch_sz; ++i)
            insert_at(U[processed + i], h_cur[i]);

        processed += cur_batch_sz;
    }

    U.clear();
}



#endif
