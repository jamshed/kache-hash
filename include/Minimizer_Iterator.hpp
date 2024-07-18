
#ifndef MINIMIZER_ITERATOR_HPP
#define MINIMIZER_ITERATOR_HPP



#include "DNA.hpp"
#include "DNA_Utility.hpp"
#include "xxHash/xxhash.h"

#include <cctype>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <cassert>


// A class to iterate over the minimizers of the constituent `k`-mers of a given
// sequence, computing the minimizer for each `k`-mer in amortized `O(1)` time.
// Canonical minimizers are computed, i.e. both strand-forms of the sequence is
// considered.
template <uint16_t k>
class Min_Iterator
{
private:

    const uint16_t l;   // Minimizer size.
    const uint16_t lmer_byte_c; // Number of bytes in an `l`-mer.
    const int16_t w;    // `l`-mer window width in a `k`-mer, i.e. number of `l`-mers in a `k`-mer.

    uint64_t last_lmer; // The last l-mer processed from the sequence.
    uint64_t last_lmer_bar; // Reverse complement of the last l-mer processed.
    const uint64_t clear_MSN_mask;  // Bitmask to clear the most-significant nucleotide bits of l-mers.

    // The hashes of the l-mers in a k-mer window. The l-mers sit at
    // `H[0..k - l]`, and is conceptually broken into two sub-windows: a prefix
    // and a suffix one.
    uint64_t H[k];

    uint64_t M_pre[k];  // Suffix-mins of the prefix window.
    uint64_t M_suf; // Min of the suffix window.

    int16_t P_pre[k];   // Positions (absolute) of the suffix-mins at the prefix window.
    int16_t P_suf;  // Position (absolute) of the min at the suffix window.

    std::size_t pivot;  // Pivot-index, inclusive to the prefix window, breaking the window into the two sub-windows.


    // Returns the minimum of `x` and `y`.
    static constexpr auto min = [](const uint64_t x, const uint64_t y){ return x < y ? x : y; };

    // Returns the 64-bit equivalent of `val`.
    template <typename T_> static constexpr uint64_t as_u64(const T_ val) { return static_cast<uint64_t>(val); }

    static constexpr auto cmov = [](const bool cond, const int16_t x, const int16_t y) { return cond ? x : y; };

    // Returns the 64-bit hash value of the canonically equivalent l-mer pair
    // `(lmer, lmer_bar)`.
    uint64_t lmer_hash(uint64_t lmer, uint64_t lmer_bar);

    // Moves the suffix window to the prefix window, and empties it—resetting
    // them to their default condition. Also recomputes the aggregate minimum
    // results of the windows.
    void reset_windows();


public:

    // Constructs a minimizer iterator to iterate over `l`-minimizers of
    // the `k`-mers of a given sequence in a streaming manner.
    Min_Iterator(uint16_t l);

    Min_Iterator(const Min_Iterator&) = delete;
    Min_Iterator(Min_Iterator&&) = delete;
    Min_Iterator& operator=(const Min_Iterator&) = delete;
    Min_Iterator& operator=(Min_Iterator&&) = delete;

    // Resets the iterator to the sequence `seq`. The iterator sits at the
    // first `k`-mer after the construction.
    void reset(const char* seq);

    // Returns the minimizer's 64-bit hash value.
    uint64_t hash() const;

    // Moves the iterator to a next k-mer extending it with the character `ch`.
    void advance(char ch);

    // Returns the minimizer's position in the `k`-mer in the forward-strand.
    // If the minimizer actually comes from the reverse-strand, then the l-mer
    // at this index must be taken in the RC-orientation to get the minimizer.
    // If there are multiple candidates, then returns some value `< 0`.
    // TODO: rename.
    int16_t hash_pos() const;
};


template <uint16_t k>
inline Min_Iterator<k>::Min_Iterator(const uint16_t l):
      l(l)
    , lmer_byte_c((2 * l + 7) / 8)
    , w(k - l + 1)
    , clear_MSN_mask(~(as_u64(0b11) << (2 * (l - 1))))
{}


template <uint16_t k>
inline uint64_t Min_Iterator<k>::lmer_hash(const uint64_t lmer, const uint64_t lmer_bar)
{
    constexpr uint64_t min_seed = 0;
    const auto h = [&](const uint64_t lmer)
        {
            // TODO: does making the size of the key unknown at compile-time non-trivially affect runtime?
            // Can use template-metaprogramming to solve this.
            return XXH3_64bits_withSeed(&lmer, lmer_byte_c, min_seed);
            // return XXH3_64bits_withSeed(&lmer, sizeof(lmer), min_seed);
        };

    return min(h(lmer), h(lmer_bar));
}


template <uint16_t k>
inline void Min_Iterator<k>::reset(const char* const seq)
{
    last_lmer = 0;
    last_lmer_bar = 0;

    for(std::size_t idx = 0; idx < l; ++idx)
    {
        assert(DNA_Utility::is_DNA_base(seq[idx]));

        const auto b = DNA_Utility::map_base(seq[idx]);
        assert(b != DNA::N);
        last_lmer |= (as_u64(b) << (2 * (l - 1 - idx)));
        last_lmer_bar |= (as_u64(DNA_Utility::complement(b)) << (2 * idx));
    }

    pivot = 0;
    H[pivot] = lmer_hash(last_lmer, last_lmer_bar);

    for(std::size_t idx = l; idx < k; ++idx)
    {
        assert(DNA_Utility::is_DNA_base(seq[idx]));

        const DNA::Base b = DNA_Utility::map_base(seq[idx]);
        assert(b != DNA::N);
        last_lmer = ((last_lmer & clear_MSN_mask) << 2) | b;
        last_lmer_bar = (last_lmer_bar >> 2) | (as_u64(DNA_Utility::complement(b)) << (2 * (l - 1)));

        H[++pivot] = lmer_hash(last_lmer, last_lmer_bar);
    }


    reset_windows();
}


template <uint16_t k>
inline void Min_Iterator<k>::reset_windows()
{
    assert(pivot == k - l);

    // TODO: try vectorization, if not already done by the compiler.

    M_pre[w] = std::numeric_limits<uint64_t>::max();
    P_pre[w] = -w;
    for(int16_t i = w - 1; i >= 0; --i)
        M_pre[i] = min(M_pre[i + 1], H[i]),
        P_pre[i] = cmov(M_pre[i + 1] < H[i], P_pre[i + 1], cmov(M_pre[i + 1] > H[i], i, -w));

    M_suf = std::numeric_limits<uint64_t>::max();
    P_suf = -w;
    pivot = 0;
}


template <uint16_t k>
inline uint64_t Min_Iterator<k>::hash() const
{
    return min(M_pre[pivot], M_suf);
}


template <uint16_t k>
inline int16_t Min_Iterator<k>::hash_pos() const
{
    return cmov(M_pre[pivot] < M_suf,
                    P_pre[pivot] - pivot,
                    cmov(M_pre[pivot] > M_suf,
                        w - pivot + P_suf, -w));
}


template <uint16_t k>
inline void Min_Iterator<k>::advance(const char ch)
{
    assert(DNA_Utility::is_DNA_base(ch));

    const DNA::Base b = DNA_Utility::map_base(ch);
    assert(b != DNA::N);
    last_lmer = ((last_lmer & clear_MSN_mask) << 2) | b;
    last_lmer_bar = (last_lmer_bar >> 2) | (as_u64(DNA_Utility::complement(b)) << (2 * (l - 1)));

    H[pivot] = lmer_hash(last_lmer, last_lmer_bar);
    P_suf = cmov(M_suf < H[pivot], P_suf, cmov(M_suf > H[pivot], pivot, -w));
    M_suf = min(M_suf, H[pivot]);

    if(pivot < static_cast<std::size_t>(k - l))
        pivot++;
    else
        reset_windows();
}



#endif
