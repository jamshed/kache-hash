
#ifndef STREAMING_KMER_HASHSET_HPP
#define STREAMING_KMER_HASHSET_HPP



#include "Kmer.hpp"
#include "Minimizer_Iterator.hpp"
#include "Directed_Vertex.hpp"
#include "utility.hpp"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <cmath>


// A fixed-sized hashset for `k`-mers.
template <uint16_t k>
class Streaming_Kmer_Hashset
{
private:

    static constexpr double lf_default = 0.75;  // Default maximum load-factor.

    const std::size_t capacity_;    // True capacity of the set; adjusted to be a power of 2.
    const std::size_t wrapper_mask; // Bitmask to wrap indexing into the table of the set.

    Kmer<k>* const K;   // Key collection.
    uint8_t* const T;   // Timestamp of the table version where each entry belongs to.

    std::size_t sz; // Size of the hashset.

    uint8_t cur_stamp;  // Current timestamp, or the version of the table. It need not be just 1 byte.

    static constexpr uint16_t l = 12;
    Min_Iterator<k> min_it;


    // Maps the hash value `h` to an index into the table.
    auto hash_to_idx(const std::size_t h) const { return h & wrapper_mask; }

    // Returns the next (wrapped) index to `i`.
    auto next_idx(const std::size_t idx) const { return hash_to_idx(idx + 1); }

    // Attemps to insert the `k`-mer `key` at the index `idx`. Returns `true`
    // iff the insertion was successful, i.e. when the `k`-mer was not present
    // at the index.
    bool try_insert_at(const Kmer<k>& key, std::size_t idx);

    // Attemps to insert the `k`-mer `key` onward from the index `idx`. Returns
    // `true` iff the insertion was successful, i.e. when the `k`-mer was not
    // present in the set.
    bool try_insert_from(const Kmer<k>& key, std::size_t idx);

    // Returns `true` iff the current key at index `idx` matches `key`.
    bool probe_at(const Kmer<k>& key, std::size_t idx);


public:

    // Constructs a hashset to support upto `max_n` `k`-mers, with a maximum
    // load-factor of `lf`.
    Streaming_Kmer_Hashset(std::size_t max_n, double lf = lf_default);

    Streaming_Kmer_Hashset(const Streaming_Kmer_Hashset&) = delete;
    Streaming_Kmer_Hashset(Streaming_Kmer_Hashset&&) = delete;
    Streaming_Kmer_Hashset& operator=(const Streaming_Kmer_Hashset&) = delete;
    Streaming_Kmer_Hashset& operator=(Streaming_Kmer_Hashset&&) = delete;

    ~Streaming_Kmer_Hashset() { deallocate(K), deallocate(T); }

    // Returns the size of the hashset.
    auto size() const { return sz; }

    // Returns the true capacity of the hashset.
    auto capacity() const { return capacity_; }

    // Clears the hashset.
    void clear();

    // Inserts the `k`-mers from the sequence `seq` of length `len` to the
    // hashset.
    void insert(const char* seq, std::size_t len);
};


template <uint16_t k>
inline Streaming_Kmer_Hashset<k>::Streaming_Kmer_Hashset(const std::size_t max_n, const double lf):
      capacity_(ceil_pow_2(static_cast<std::size_t>(std::ceil(max_n / lf))))
    , wrapper_mask(capacity_ - 1)
    , K(allocate<Kmer<k>>(capacity_))
    , T(allocate<uint8_t>(capacity_))
    , sz(0)
    , cur_stamp(std::numeric_limits<uint8_t>::max())
    , min_it(l)
{
    assert(capacity_ > 1);

    clear();
}


template <uint16_t k>
inline void Streaming_Kmer_Hashset<k>::clear()
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
inline void Streaming_Kmer_Hashset<k>::insert(const char* const seq, const std::size_t len)
{
    assert(len >= k);

    Directed_Vertex<k> v(seq);  // Current vertex in the scan over `seq`.
    std::size_t kmer_idx = 0;   // Index of the current `k`-mer.

    min_it.reset(seq);
    while(true)
    {
        const auto h = min_it.hash();
        const auto p = min_it.hash_pos();
        const auto key = v.canonical();

        if(p < 0)   // Ambiguous minimizer-instance.
            // `k`-mer gets hashed independently.
            try_insert_from(key, hash_to_idx(key.to_u64()));
        else    // Unambiguous minimizer-instance.
        {
            // Hash the minimizer and dock `k`-mer accordingly.
            const auto d1 = hash_to_idx(h + p);
            const auto d2 = hash_to_idx(h + (k - (p + l)));
            if(!probe_at(key, d1) && !probe_at(key, d2))
                if(!try_insert_at(key, d1) && !try_insert_at(key, d2))
                    try_insert_from(key, hash_to_idx(key.to_u64()));
        }

        if(kmer_idx + k == len)
            break;

        min_it.advance(seq[kmer_idx + k]);
        v.roll_forward(DNA_Utility::map_base(seq[kmer_idx + k]));
        kmer_idx++;
    }
}


template <uint16_t k>
inline bool Streaming_Kmer_Hashset<k>::try_insert_at(const Kmer<k>& key, const std::size_t idx)
{
    if(T[idx] != cur_stamp) // Empty slot as it contains an entry from a previous version of the table.
    {
        assert(size() < capacity());

        T[idx] = cur_stamp;
        K[idx] = key;
        sz++;

        return true;
    }

    return false;
}


template <uint16_t k>
inline bool Streaming_Kmer_Hashset<k>::try_insert_from(const Kmer<k>& key, const std::size_t idx)
{
    for(std::size_t i = idx; ; i = next_idx(i))
        if(try_insert_at(key, i))
            return true;
        else if(K[i] == key)
            return false;

    return false;
}


template <uint16_t k>
inline bool Streaming_Kmer_Hashset<k>::probe_at(const Kmer<k>& key, const std::size_t idx)
{
    return T[idx] == cur_stamp && K[idx] == key;
}



#endif
