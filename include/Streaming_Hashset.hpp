
#ifndef STREAMING_HASHSET_HPP
#define STREAMING_HASHSET_HPP



#include "Kmer.hpp"
#include "DNA_Utility.hpp"
#include "Directed_Vertex.hpp"
#include "Minimizer_Iterator.hpp"
#include "utility.hpp"
#include "unordered_dense/unordered_dense.h"

#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <cstddef>
#include <cstring>
#include <cassert>
#include <vector>


// Hash set particularly suited for streaming k-mer operations.
template <uint16_t k, uint16_t l = 19>
class Streaming_Hash_Set
{
private:

    static constexpr double lf_default = 0.375; // Default maximum load-factor.

    static constexpr std::size_t B = 8;    // Block size.

    std::size_t cap_;   // True capacity of the set in terms of blocks; adjusted to be a power of 2.
    std::size_t idx_mask;   // Mask to wrap indexing within the table.

    Kmer<k> (* K)[B];   // Key collection, blocked.
    uint8_t (* T)[B];   // Timestamp of the table version where each entry belongs to.

    std::size_t sz; // Size of the set.

    uint8_t cur_stamp;  // Current timestamp, i.e. the version of the table.

    ankerl::unordered_dense::set<Kmer<k>, Kmer_Hasher<k>> overflow; // The overflow set.

    uint64_t sum_len;   // Total sequence length provided.
    uint64_t sup_kmer_c;    // Number of super k-mers in the provided sequences to the set.

    ankerl::unordered_dense::set<uint64_t, Kmer_Hasher<k>> mins; // Set of all minimizer; used for developer stats.
    ankerl::unordered_dense::set<uint64_t, Kmer_Hasher<k>> mins_ov;  // Set of overflowing minimizers; used for developer stats.
    std::vector<uint64_t> M;


public:

    // Constructs a hash set to support upto `max_sz` `k`-mers, with a maximum
    // load factor of `lf`.
    Streaming_Hash_Set(std::size_t max_sz, double lf = lf_default);

    Streaming_Hash_Set(const Streaming_Hash_Set&) = delete;
    Streaming_Hash_Set(Streaming_Hash_Set&&) = delete;
    Streaming_Hash_Set& operator=(const Streaming_Hash_Set&) = delete;
    Streaming_Hash_Set& operator=(Streaming_Hash_Set&&) = delete;

    ~Streaming_Hash_Set() { deallocate(K), deallocate(T); }

    // Returns the capacity of the hash set.
    auto capacity() const { return cap_ * B; }

    // Returns the size of the hash set.
    auto size() const { return sz + overflow.size(); }

    // Clears the hash set.
    void clear();

    // Inserts `k`-mers from the `len`-length sequence `s` into the set.
    template <bool track_stats_ = false>
    void insert(const char* s, std::size_t len);

    // Prints developer statistics of the hash set.
    void print_stats();

    // TODO: flesh out API.
};


template <uint16_t k, uint16_t l>
inline Streaming_Hash_Set<k, l>::Streaming_Hash_Set(const std::size_t max_sz, const double lf):
      cap_(ceil_pow_2(static_cast<std::size_t>((static_cast<double>(max_sz) / B) / lf)))
    , idx_mask(cap_ - 1)
    , K(allocate<Kmer<k>[B]>(cap_))
    , T(allocate<uint8_t[B]>(cap_))
    , cur_stamp(std::numeric_limits<uint8_t>::max())
{
    assert(capacity() > 1);

    clear();
}


template <uint16_t k, uint16_t l>
inline void Streaming_Hash_Set<k, l>::clear()
{
    sz = 0;
    sum_len = 0;
    sup_kmer_c = 0;
    overflow.clear();

    cur_stamp = (cur_stamp + 1u) & std::numeric_limits<decltype(cur_stamp)>::max();
    if(cur_stamp == 0)
    {
        std::memset(reinterpret_cast<char*>(T), 0, capacity() * sizeof(cur_stamp));
        cur_stamp = 1;
    }

    M.resize(cap_, 0);  // TODO: only perform when tracking with dev-stats.
}


template <uint16_t k, uint16_t l>
template <bool track_stats_>
inline void Streaming_Hash_Set<k, l>::insert(const char* const s, const std::size_t len)
{
    sum_len += len;
    if(len < k)
        return;

    Directed_Vertex<k> v(s);
    std::size_t kmer_idx = 0;

    Min_Iterator<k> min_it(l);
    min_it.reset(s);
    uint64_t last_min = 0;  // Undercounting 1 super k-mer from `s` with probability 1/2^64.
    while(true)
    {
        const auto h = min_it.hash();
        const auto key = v.canonical();

        if constexpr(track_stats_)
        {
            sup_kmer_c += (h != last_min);
            last_min = h;

            mins.insert(h);
        }


        const auto i = h & idx_mask;
        std::size_t j = 0;

        if constexpr(track_stats_)
        {
            if(M[i] == 0)
                M[i] = h;
        }

        for(; j < B; ++j)
            if(T[i][j] != cur_stamp)
            {
                T[i][j] = cur_stamp;
                K[i][j] = key;
                sz++;

                break;
            }
            else if(K[i][j] == key)
                break;

        if(j == B)
        {
            if constexpr(track_stats_)
                mins_ov.insert(h);

            overflow.insert(key);
        }

        if(kmer_idx + k == len)
            break;

        min_it.advance(s[kmer_idx + k]);
        v.roll_forward(DNA_Utility::map_base(s[kmer_idx + k]));
        kmer_idx++;
    }
}


template <uint16_t k, uint16_t l>
inline void Streaming_Hash_Set<k, l>::print_stats()
{
    std::cerr   << "Total sequence length:  " << sum_len << ".\n";
    std::cerr   << "Number of super k-mers: " << sup_kmer_c << ".\n";
    std::cerr   << "========================================\n";

    uint64_t hash_conflict = 0;
    for(const auto& m : mins)
        if(M[m & idx_mask] != m)
            hash_conflict++;

    std::cerr   << "Minimizer set size:          " << mins.size() << ".\n";
    std::cerr   << "Overflow minimizer set size: " << mins_ov.size() << " ("
                    << (mins_ov.size() * 100.0 / mins.size()) << "%).\n";
    std::cerr   << "Number of hash-conflicts for minimizers: " << hash_conflict
                    << " (" << (hash_conflict * 100.0 / mins.size()) << "%).\n";
    std::cerr   << "========================================\n";

    std::cerr   << "Complete set size: " << size() << ".\n";
    std::cerr   << "Overflow set size: " << overflow.size() << " ("
                    << (overflow.size() * 100.0 / size()) << "%).\n";
    std::cerr   << "========================================\n";

    uint64_t full = 0, empty = 0;
    for(std::size_t i = 0; i < cap_; ++i)
    {
        std::size_t j = 0;
        for(; j < B; ++j)
            if(T[i][j] != cur_stamp)
                break;

        full += (j == B);
        empty += (j == 0);
    }

    std::cerr   << "#Blocks:       " << cap_ << ".\n";
    std::cerr   << "#Full-blocks:  " << full << "(" << (full * 100.0 / cap_) << "%).\n";
    std::cerr   << "#Empty-blocks: " << empty << "(" << (empty * 100.0 / cap_) << "%).\n";
    std::cerr   << "========================================\n";
}



#endif
