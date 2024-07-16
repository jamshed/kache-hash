

#include "Kmer.hpp"

#include <cstdint>
#include <cstddef>


template <typename T_seq_, uint16_t k, bool is_canonical_> class Minimizer_Iterator;


// =============================================================================
// A class containing various utility methods for k-mer minimizers.
class Minimizer_Utility
{
    typedef uint64_t minimizer_t;

public:

    // Returns the hash value of the l-mer `lmer`. The seed-value `seed` is used
    // in hashing.
    static uint64_t hash(minimizer_t lmer, uint64_t seed = 0);

    // Extracts the minimizer of the k-mer `kmer` into `min`, and its index into
    // `idx`.
    template <uint16_t k>
    static void get_minimizer(const Kmer<k>& kmer, const uint16_t l, minimizer_t& min, std::size_t& idx);

    // Extracts the minimizer of the k-mer `kmer` into `min`.
    template <uint16_t k>
    static minimizer_t canonical_minimizer(const Kmer<k>& kmer, const uint16_t l);
};


// A class to pack information regarding l-mers, to aid in computing
// l-minimizers of k-mers over sequences by `Minimizer_Iterator`.
class Lmer_Tuple
{
    friend class Minimizer_Utility;
    template <typename T_seq_, uint16_t k, bool is_canonical_> friend class Minimizer_Iterator;

    typedef uint64_t minimizer_t;

private:

    minimizer_t lmer;   // The l-mer.
    std::size_t index;  // Its index in the underlying sequence.
    uint64_t hash;      // Its hash value, determining the l-mer ordering.


public:

    Lmer_Tuple()
    {};

    // Constructs a tuple for an l-mer `lmer`, positioned at index `index` of
    // the underlying sequence, and having a hash value `hash`.
    Lmer_Tuple(minimizer_t lmer, std::size_t index, uint64_t hash);

    // Returns `true` iff this l-mer tuple is to be ordered as lesser to the
    // tuple `rhs`. The primary ordering is based on the tuples' hashes. If
    // equal, then the ordering is based on their literal form. If also equal,
    // then the tuple occurring earlier in the underlying sequence is lesser.
    bool operator<(const Lmer_Tuple& rhs) const;
};


inline Lmer_Tuple::Lmer_Tuple(const minimizer_t lmer, const std::size_t index, const uint64_t hash):
    lmer(lmer),
    index(index),
    hash(hash)
{}


inline bool Lmer_Tuple::operator<(const Lmer_Tuple& rhs) const
{
    if(hash != rhs.hash)
        return hash < rhs.hash;

    if(lmer != rhs.lmer)
        return lmer < rhs.lmer;

    return index < rhs.index;
}


inline uint64_t Minimizer_Utility::hash(const minimizer_t lmer, const uint64_t seed)
{
#ifdef CF_DEVELOP_MODE
    return lmer;    // TODO: add as debug option for developer.
#endif

    return XXH3_64bits_withSeed(&lmer, sizeof(lmer), seed);
    // TODO: the entropy of the minimizer-hash is getting diluted by considering all the leading 0 bits (64 - 2l bits, specifically.)
    // return wyhash(&lmer, (2 * l + 7) / 8, seed, wy_salt);
    // return wyhash(&lmer, sizeof(lmer), seed, wy_salt);
}


template <uint16_t K>
inline void Minimizer_Utility::get_minimizer(const Kmer<K>& kmer, const uint16_t l, minimizer_t& min, std::size_t& idx)
{
    const uint64_t* const kmer_data = kmer.data();

    const minimizer_t last_lmer = kmer_data[0] & ((uint64_t(0b1) << (2 * l)) - 1);  // The last l-mer in `kmer`.
    Lmer_Tuple curr_lmer(last_lmer, K - l, hash(last_lmer));
    Lmer_Tuple min_lmer(curr_lmer);

    uint64_t base;
    for(std::size_t i = l; i < K; ++i)  // Backward scan of l-mers in `kmer`. Yes, backwards.
    {
        base = (kmer_data[i >> 5] >> (2 * (i & 31))) & uint64_t(0b11);
        curr_lmer.lmer = (curr_lmer.lmer >> 2) | (base << (2 * (l - 1)));
        curr_lmer.index--;
        curr_lmer.hash = hash(curr_lmer.lmer);

        if(curr_lmer < min_lmer)
            min_lmer = curr_lmer;
    }

    min = min_lmer.lmer, idx = min_lmer.index;
}


template <uint16_t k>
inline uint64_t Minimizer_Utility::canonical_minimizer(const Kmer<k>& kmer, const uint16_t l)
{
    minimizer_t min_f, min_r;
    std::size_t min_idx;

    get_minimizer(kmer, l, min_f, min_idx);
    get_minimizer(kmer.reverse_complement(), l, min_r, min_idx);

    const auto h_f = hash(min_f), h_r = hash(min_r);
    return h_f != h_r ? (h_f < h_r ? min_f : min_r) : std::min(min_f, min_r);
}



#include "utility.hpp"

#include <cstddef>


// Fixed-capacity deque of capacity at least `req_cap_`. UB when size exceeds
// capacity.
template <typename T_, std::size_t req_cap>
class Fixed_Cap_Deque
{
private:

    static constexpr std::size_t cap_ = ceil_pow_2(req_cap);    // Maximum capacity.
    static constexpr auto wrap_mask = cap_ - 1; // Mask to wrap indices around the deque.

    std::size_t front_; // Current index to the front.
    std::size_t back_;  // Current index to the back.
    T_ arr[cap_];   // Underlying memory pool.

    void grow_back() { back_ = (back_ + 1) & wrap_mask; }

    void grow_front() { front_ = (front_ - 1) & wrap_mask; }    // Under(over)-flow is defined for unsigned.

    void shrink_back() { back_ = (back_ - 1) & wrap_mask; }

    void shrink_front() { front_ = (front_ + 1) & wrap_mask; }

public:

    Fixed_Cap_Deque():
          front_(0)
        , back_(0)
    {}

    bool empty() const { return front_ == back_; }

    const T_& front() const { return arr[front_]; }

    const T_& back() const { return arr[(back_ - 1) & wrap_mask]; }

    void clear() { front_ = back_; }

    void push_back(const T_& val) { arr[back_] = val; grow_back(); }

    void push_front(const T_& val) { grow_front(); arr[front_] = val; }

    template <typename... Args>
    void emplace_back(Args&&... args) { new (arr + back_) T_(args...); grow_back(); }

    template <typename... Args>
    void emplace_front(Args&&... args) { grow_front(); new (arr + front_) T_(args...); }

    void pop_back() { shrink_back(); }

    void pop_front() { shrink_front(); }
};



#include "DNA_Utility.hpp"
#include "Kmer.hpp"
#include "xxHash/xxhash.h"

#include <cstdint>
#include <cassert>


// =============================================================================
// A class to iterate over the minimizers of the constituent `k`-mers of a given
// sequence of type `T_seq_`, computing the minimizer for each `k`-mer in
// amortized O(1) time. If `is_canonical_` is `true`, then canonical minimizers
// are computed, i.e. both strand-forms of the sequence is considered.
template <typename T_seq_, uint16_t k, bool is_canonical_ = false>
class Minimizer_Iterator
{
    typedef uint64_t minimizer_t;

private:

    T_seq_ seq; // The sequence on which to iterate over.
    std::size_t seq_len;    // Length of the given sequence.

    const uint16_t l;   // Size of the minimizers.

    const uint64_t seed;    // Seed for hashing l-mers.

    minimizer_t last_lmer;  // The last l-mer processed.
    minimizer_t last_lmer_bar;  // Reverse complement of the last l-mer processed.
    std::size_t last_lmer_idx;  // Index into the sequence of the last l-mer processed.
    const minimizer_t clear_MSN_mask;   // Bitmask to clear the most-significant nucleotide bits of l-mers.

    typedef Fixed_Cap_Deque<Lmer_Tuple, 2 * k> deque_t;

    // Collection of l-mers that have already been seen, and cannot be ruled out
    // yet as candidate minimizers for the k-mers yet to be seen fully. `dq_f`
    // is for the forward-strand and `dq_r` is for the reverse-strand.
    deque_t dq_f, dq_r;

public:

    // Constructs a minimizer iterator to iterate over `l`-minimizers of
    // the `k`-mers of a given sequence in a streaming manner. The seed-value
    // `seed` is used in hashing the `l`-mers.
    Minimizer_Iterator(uint16_t l, uint64_t seed = 0);

    // Constructs a minimizer iterator to iterate over `l`-minimizers of
    // the `k`-mers of the sequence `seq`, of length `seq_len`, in a streaming
    // manner. The iterator sits at the first k-mer after the construction.
    // The seed-value `seed` is used in hashing the `l`-mers.
    Minimizer_Iterator(T_seq_ seq, std::size_t seq_len, uint16_t l, uint64_t seed = 0);

    // Resets the iterator to the sequence `seq` of length `seq_len`. The
    // iterator sits at the first k-mer after the construction.
    void reset(T_seq_ seq, std::size_t seq_len);

    // Moves the iterator to a next k-mer by extending it with the character
    // `ch`.
    void advance(char ch);

    // Moves the iterator to the next k-mer in the sequence. Returns `true` iff
    // the current k-mer is not the last k-mer in the sequence.
    bool operator++();

    // Puts the minimizer of the current k-mer into `minimizer`, and stores its
    // index in the sequence at `index`. The position information over the
    // sequence is maintained internally.
    void value_at(minimizer_t& minimizer, std::size_t& index) const;

    // Puts the minimizer of the current k-mer into `minimizer`, stores its
    // index in the sequence at `index`, and stores its 64-bit hash at `h`. The
    // position information over the sequence is maintained internally.
    void value_at(minimizer_t& minimizer, std::size_t& index, uint64_t& h) const;

    // Returns the minimizer's 64-bit hash value.
    uint64_t hash() const;

    // Computes the `l`-minimizer of the `k`-length sequence `kmer` at `min`,
    // its hash at `h`, and index in `kmer` ar `idx`. The seed value `seed` is
    // used in hashing the `l`-mers.
    static void minimizer(T_seq_ kmer, uint16_t l, uint64_t seed, minimizer_t& min, uint64_t& h, std::size_t& idx);
};


template <typename T_seq_,  uint16_t k, bool is_canonical_>
inline Minimizer_Iterator<T_seq_, k, is_canonical_>::Minimizer_Iterator(const uint16_t l, const uint64_t seed):
      l(l)
    , seed(seed)
    , clear_MSN_mask(~(static_cast<minimizer_t>(0b11) << (2 * (l - 1))))
{
    assert(l <= k);
}


template <typename T_seq_, uint16_t k, bool is_canonical_>
inline Minimizer_Iterator<T_seq_, k, is_canonical_>::Minimizer_Iterator(T_seq_ const seq, const std::size_t seq_len, const uint16_t l, const uint64_t seed):
    Minimizer_Iterator(l, seed)
{
    reset(seq, seq_len);
}


template <typename T_seq_, uint16_t k, bool is_canonical_>
inline void Minimizer_Iterator<T_seq_, k, is_canonical_>::reset(T_seq_ const seq, const std::size_t seq_len)
{
    this->seq = seq;
    this->seq_len = seq_len;
    dq_f.clear(), dq_r.clear();

    last_lmer = 0;
    last_lmer_bar = 0;
    last_lmer_idx = 0;

    DNA::Base base;
    for(std::size_t idx = 0; idx < l; ++idx)
    {
        base = DNA_Utility::map_base(seq[idx]);
        assert(base != DNA::Base::N);
        last_lmer |= (static_cast<minimizer_t>(base) << (2 * (l - 1 - idx)));

        if constexpr(is_canonical_)
            last_lmer_bar |= (static_cast<minimizer_t>(DNA_Utility::complement(base)) << (2 * idx));
    }

    dq_f.emplace_back(last_lmer, last_lmer_idx, Minimizer_Utility::hash(last_lmer, seed));
    if constexpr(is_canonical_)
        dq_r.emplace_back(last_lmer_bar, last_lmer_idx, Minimizer_Utility::hash(last_lmer_bar, seed));

    while(last_lmer_idx + (l - 1) < static_cast<std::size_t>(k - 1))
        operator++();
}


template <typename T_seq_, uint16_t k, bool is_canonical_>
inline void Minimizer_Iterator<T_seq_, k, is_canonical_>::advance(const char ch)
{
    assert(DNA_Utility::is_DNA_base(ch));

    last_lmer_idx++;
    const DNA::Base base = DNA_Utility::map_base(ch);
    last_lmer = ((last_lmer & clear_MSN_mask) << 2) | base;
    if constexpr(is_canonical_)
        last_lmer_bar = (last_lmer_bar >> 2) | (static_cast<minimizer_t>(DNA_Utility::complement(base)) << (2 * (l - 1)));


    // End of the current k-mer, e = last_lmer_idx + l - 1; So start of the current k-mer: e - (k - 1).
    const std::size_t curr_kmer_idx = (last_lmer_idx + l - 1 >= k - 1lu ? last_lmer_idx + l - 1 - (k - 1) : 0);
    if(dq_f.front().index < curr_kmer_idx)    // This candidate l-mer falls out of the current k-mer.
        dq_f.pop_front();

    if constexpr(is_canonical_)
        if(dq_r.front().index < curr_kmer_idx)  // This candidate l-mer falls out of the current k-mer.
            dq_r.pop_front();


    const auto fix_dq =
        [](deque_t& dq, const Lmer_Tuple& lmer_tuple)
        {
            while(!dq.empty())
                if(dq.back() < lmer_tuple)
                    break;
                else
                    dq.pop_back();

            dq.push_back(lmer_tuple);
        };


    fix_dq(dq_f, Lmer_Tuple(last_lmer, last_lmer_idx, Minimizer_Utility::hash(last_lmer, seed)));
    if constexpr(is_canonical_)
        fix_dq(dq_r, Lmer_Tuple(last_lmer_bar, last_lmer_idx, Minimizer_Utility::hash(last_lmer_bar, seed)));
}


template <typename T_seq_, uint16_t k, bool is_canonical_>
inline bool Minimizer_Iterator<T_seq_, k, is_canonical_>::operator++()
{
    if(last_lmer_idx + l == seq_len)
        return false;

    advance(seq[last_lmer_idx + l]);
    return true;
}


template <typename T_seq_, uint16_t k, bool is_canonical_>
inline void Minimizer_Iterator<T_seq_, k, is_canonical_>::value_at(minimizer_t& minimizer, std::size_t& index) const
{
    if constexpr(is_canonical_)
    {
        if(dq_f.front() < dq_r.front())
            minimizer = dq_f.front().lmer, index = dq_f.front().index;
        else
            minimizer = dq_r.front().lmer, index = dq_r.front().index;
    }
    else
        minimizer = dq_f.front().lmer, index = dq_f.front().index;
}


template <typename T_seq_, uint16_t k, bool is_canonical_>
inline void Minimizer_Iterator<T_seq_, k, is_canonical_>::value_at(minimizer_t& minimizer, std::size_t& index, uint64_t& h) const
{
    if constexpr(is_canonical_)
    {
        if(dq_f.front() < dq_r.front())
            minimizer = dq_f.front().lmer, index = dq_f.front().index, h = dq_f.front().hash;
        else
            minimizer = dq_r.front().lmer, index = dq_r.front().index, h = dq_r.front().hash;
    }
    else
        minimizer = dq_f.front().lmer, index = dq_f.front().index, h = dq_f.front().hash;
}


template <typename T_seq_, uint16_t k, bool is_canonical_>
inline uint64_t Minimizer_Iterator<T_seq_, k, is_canonical_>::hash() const
{
    if constexpr(is_canonical_)
        return dq_f.front() < dq_r.front() ? dq_f.front().hash : dq_r.front().hash;
    else
        return dq_f.front().hash;
}


template <typename T_seq_, uint16_t k, bool is_canonical_>
inline void Minimizer_Iterator<T_seq_, k, is_canonical_>::minimizer(T_seq_ kmer, uint16_t l, uint64_t seed, minimizer_t& min, uint64_t& h, std::size_t& idx)
{
    Minimizer_Iterator(kmer, k, l, seed).value_at(min, idx);
    h = Minimizer_Utility::hash(min);
}
