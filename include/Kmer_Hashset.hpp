
#ifndef KMER_HASHSET_HPP
#define KMER_HASHSET_HPP



#include "Kmer.hpp"
#include "utility.hpp"

#include <cstdint>
#include <cstddef>
#include <limits>
#include <cstring>
#include <cmath>
#include <cassert>


// A fixed-sized hashset for k-mers.
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


    // Maps the hash value `h` to an index into the table.
    auto hash_to_idx(const std::size_t h) const { return h & wrapper_mask; }

    // Returns the next (wrapped) index to `i`.
    auto next_idx(const std::size_t idx) const { return hash_to_idx(idx + 1); }


public:

    // Constructs a hashset to support upto `max_n` k-mers, with a maximum
    // load-factor of `lf`.
    Kmer_Hashset(std::size_t max_n, double lf = lf_default);

    Kmer_Hashset(const Kmer_Hashset&) = delete;
    Kmer_Hashset(Kmer_Hashset&&);
    Kmer_Hashset& operator=(const Kmer_Hashset&) = delete;
    Kmer_Hashset& operator=(Kmer_Hashset&&) = delete;

    ~Kmer_Hashset() { deallocate(T); }

    // Returns the size of the hashset.
    auto size() const { return sz; }

    // Returns the true capacity of the hashset.
    auto capacity() const { return capacity_; }

    // Clears the hashset.
    void clear();

    // Inserts the k-mer `key` to the hashset. Returns `true` iff the insertion
    // was successful, i.e. when the k-mer was not present in the set.
    bool insert(const Kmer<k>& key);
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
    const auto h = hash_to_idx(key.to_u64());
    for(std::size_t i = h; ; i = next_idx(i))
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



#endif
