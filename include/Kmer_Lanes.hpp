
#ifndef KACHE_HASH_KMER_LANES_HPP
#define KACHE_HASH_KMER_LANES_HPP



#include "Kmer_Batch.hpp"

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <bit>
#include <algorithm>
#include <utility>
#include <cassert>

#include <immintrin.h>


namespace kache_hash
{


// Describes the `k`-mers of up to `L` stretches of ACGT bases at once, one per
// SIMD lane, in lock-step chunks of `T` bases: what `Kmer_Batch` gives for each
// stretch, i.e. the canonical forms, checksums, and `l`-minimizer hashes.
template <uint16_t k, uint16_t l>
class Kmer_Lanes
{
    static_assert(k <= 32, "Lane-parallel k-mer descriptions need k-mers in a machine word.");
    static_assert(l <= k && l <= 32);

public:

    static constexpr std::size_t L = 4; // Number of lanes: 4 keeps the vectors at 256 bits, which spares the clock frequency.
    static constexpr std::size_t T = 16;    // Number of bases a chunk advances each lane by.

    typedef uint64_t vec __attribute__((vector_size(32)));  // A 64-bit value per lane.

    // Descriptions of a chunk of k-mers, a vector per base: `key[t][j]` etc.
    // are of the k-mer ending at the chunk's `t`'th base in the `j`'th lane,
    // which is the `first + t`'th k-mer of the lane's stretch; it exists iff
    // that index is in `[0, n[j])`. Stored as computed, they need no transpose.
    struct Chunk
    {
        vec key[T]; // Canonical forms.
        vec h[T];   // Minimizer hashes.
        vec c[T];   // Checksums, never 0, in the low byte; bit 8 is set iff the k-mer, as it occurs, is its canonical
                    // form, when the orientations are asked for.
        int64_t first;  // Index of the k-mers at the chunk's first base, in their stretches; negative before the first k-mers.
        std::size_t n[L];   // Number of k-mers in the lanes' stretches.
        bool any;   // Whether the chunk has any k-mer.

        // Returns the bounds `[t_b, t_e)` of the bases ending the `j`'th lane's k-mers.
        std::pair<std::size_t, std::size_t> range(const std::size_t j) const
        {
            const int64_t b = std::max<int64_t>(0, -first), e = std::min<int64_t>(T, static_cast<int64_t>(n[j]) - first);
            return {static_cast<std::size_t>(b), static_cast<std::size_t>(std::max(b, e))};
        }
    };

private:

    typedef int64_t ivec __attribute__((vector_size(32)));  // A lane-wise mask: all ones or zero.
    typedef vec vec_unaligned __attribute__((may_alias, aligned(8)));   // For loading lanes from 64-bit words.

    static constexpr uint16_t w = k - l + 1;    // Number of l-mers in a k-mer.
    static constexpr std::size_t P = T + w - 1; // Number of l-mer positions spanned by a chunk's k-mers.
    static constexpr uint64_t kmer_mask = (k == 32 ? ~uint64_t(0) : (uint64_t(1) << (2 * k)) - 1);
    static constexpr uint64_t lmer_mask = (l == 32 ? ~uint64_t(0) : (uint64_t(1) << (2 * l)) - 1);
    static constexpr uint16_t lvl = std::bit_width(w) - 1;  // Number of doubling levels of the sliding minimum.
    static constexpr std::size_t span = std::size_t(1) << lvl;
    static constexpr std::size_t off = w - span;    // A window `[i, i + w)` is the union of the spans at `i` and `i + off`.

    const char* t[L];   // The stretches.
    std::size_t len[L]; // Their lengths, in bases.
    std::size_t n[L];   // Their numbers of k-mers.
    static constexpr int64_t pad = (T - (k - 1) % T) % T;  // Number of bases read before the stretches, as A's, so that the first k-mers end at a chunk's first base: a stretch's k-mers then fill whole chunks, rather than starting with a few.

    int64_t s;  // Number of bases each lane has been advanced by, less `pad`.
    std::size_t s_end;  // Most bases any lane advances by.

    vec kf, kr;     // Forward and reverse forms of the lanes' last k-mers, or of their prefixes.
    vec lf, lr;     // The same of their last l-mers.

    alignas(64) vec hc[P];  // `hc[x]` is the hash of the canonical form of the lanes' l-mers at the chunk's `x`'th position; the first `w - 1` are the last chunk's.

    // Returns `x` broadcast to all lanes.
    static vec splat(const uint64_t x) { return vec{} + x; }

    // Returns the lane-wise selection of `a` where `c` is all ones, else `b`.
    static vec select(const ivec c, const vec a, const vec b) { return c ? a : b; }

    // Returns the lane-wise minimum of `a` and `b`.
    static vec min(const vec a, const vec b) { return a < b ? a : b; }

    // Returns the codes of the next `T` bases of the lanes, 8 bits per base.
    void load_bases(vec (&codes)[T / 8]) const;

public:

    // Resets the describer to the stretches `t[0..c)` of lengths `len[0..c)`,
    // each of at least `k` bases; the remaining lanes idle.
    void reset(const char* const* t, const std::size_t* len, std::size_t c);

    // Describes the next chunk into `out`, with the k-mers' orientations iff
    // `orient_`. Returns `false` iff all lanes are exhausted; then `out` holds
    // no k-mer.
    template <bool orient_ = false> bool next(Chunk& out);

    // Puts the minimizer hashes of the `n` canonical k-mers with the encodings
    // `x[0..n)` at `h[0..n)`, as the describers give them, from all their
    // l-mers: a k-mer per lane.
    static void minimizer_hashes(const uint64_t* x, std::size_t n, uint64_t* h);
};


template <uint16_t k, uint16_t l>
inline void Kmer_Lanes<k, l>::reset(const char* const* const t, const std::size_t* const len, const std::size_t c)
{
    assert(c <= L);
    s = -pad, s_end = 0;
    for(std::size_t j = 0; j < L; ++j)
    {
        this->t[j] = (j < c ? t[j] : nullptr);
        this->len[j] = (j < c ? len[j] : 0);
        assert(this->len[j] == 0 || this->len[j] >= k);
        n[j] = (this->len[j] >= k ? this->len[j] - (k - 1) : 0);
        s_end = std::max(s_end, this->len[j]);
    }

    kf = kr = lf = lr = splat(0);
}


template <uint16_t k, uint16_t l>
inline void Kmer_Lanes<k, l>::load_bases(vec (&codes)[T / 8]) const
{
    static_assert(T % 8 == 0);
    for(std::size_t g = 0; g < T / 8; ++g)
    {
        alignas(64) uint64_t x[L];
        for(std::size_t j = 0; j < L; ++j)
        {
            const int64_t b = s + 8 * static_cast<int64_t>(g);  // First base of the group.
            const auto n_j = static_cast<int64_t>(len[j]);
            if(b >= 0 && b + 8 <= n_j)
                std::memcpy(&x[j], t[j] + b, 8);
            else
            {
                x[j] = 0;
                const auto f = std::max<int64_t>(b, 0), e = std::min<int64_t>(b + 8, n_j);  // The group's bases within the stretch.
                if(f < e)
                    std::memcpy(reinterpret_cast<char*>(&x[j]) + (f - b), t[j] + f, e - f);
            }
        }

        const vec v = *reinterpret_cast<const vec_unaligned*>(x);
        codes[g] = ((v >> 1) ^ (v >> 2)) & splat(0x0303030303030303ull);  // The 2-bit code of each base, in its byte.
    }
}


template <uint16_t k, uint16_t l>
template <bool orient_>
inline bool Kmer_Lanes<k, l>::next(Chunk& out)
{
    out.any = false;
    if(s >= static_cast<int64_t>(s_end))
        return false;

    vec codes[T / 8];
    load_bases(codes);

    // The lanes advance together, so a chunk ends no k-mer, or completes no l-mer in a k-mer's window, in all lanes
    // at once: the first two chunks of reads do neither, or only the latter. These chunks skip the work for them.
    const bool lmers = (s + static_cast<int64_t>(T) > static_cast<int64_t>(l) - 1);    // Whether the chunk completes l-mers.
    const bool kmers = (s - static_cast<int64_t>(k - 1) + static_cast<int64_t>(T) > 0);    // Whether the chunk ends k-mers.

    // Roll the k-mers and the l-mers by the chunk's bases, hashing the l-mers. Before a stretch's first l-mer, the
    // positions hold garbage, which no k-mer's window covers.
    alignas(64) vec f[T], r[T];  // The k-mers' forms after each base.
    for(std::size_t b = 0; b < T; ++b)
    {
        const vec base = (codes[b / 8] >> (8 * (b % 8))) & splat(3);
        const vec comp = base ^ splat(3);
        kf = ((kf << 2) | base) & splat(kmer_mask);
        kr = (kr >> 2) | (comp << (2 * (k - 1)));
        lf = ((lf << 2) | base) & splat(lmer_mask);
        lr = (lr >> 2) | (comp << (2 * (l - 1)));
        f[b] = kf, r[b] = kr;

        if(lmers)
            hc[w - 1 + b] = wang_u64(min(lr, lf));
    }

    if(!kmers)
    {
        out.first = s - (k - 1);
        for(std::size_t x = 0; x + T < P; ++x)
            hc[x] = hc[x + T];
        for(std::size_t j = 0; j < L; ++j)
            out.n[j] = n[j];
        s += static_cast<int64_t>(T);
        return true;
    }

    // Sliding minima over `w` positions, by doubling spans.
    alignas(64) vec v[P];
    for(std::size_t x = 0; x < P; ++x)
        v[x] = hc[x];
    std::size_t sp = 1;
    for(uint16_t j = 0; j < lvl; ++j, sp *= 2)
        for(std::size_t x = 0; x + 2 * sp <= P; ++x)
            v[x] = min(v[x], v[x + sp]);

    // The k-mers ending at the chunk's bases: the `i`'th of them ends at base `s + i`, and its window is at
    // positions `[i, i + w)` of `hc`.
    out.first = s - (k - 1);
    for(std::size_t i = 0; i < T; ++i)
    {
        const vec h = min(v[i], v[i + off]);

        const ivec canon_fwd = (f[i] < r[i]);
        const vec key = select(canon_fwd, f[i], r[i]);
        const vec c = checksum_hash(key) >> 56; // As `kmer_checksum`.

        out.key[i] = key;
        out.h[i] = h;
        out.c[i] = select(c == splat(0), splat(1), c);
        if constexpr(orient_)
            out.c[i] |= select(key == f[i], splat(0x100), splat(0));    // A palindrome counts as forward.
    }

    // Keep the last `w - 1` positions.
    for(std::size_t x = 0; x + T < P; ++x)
        hc[x] = hc[x + T];

    for(std::size_t j = 0; j < L; ++j)
        out.n[j] = n[j];
    out.any = (out.first + static_cast<int64_t>(T) > 0);
    s += static_cast<int64_t>(T);
    return true;
}


template <uint16_t k, uint16_t l>
inline void Kmer_Lanes<k, l>::minimizer_hashes(const uint64_t* const x, const std::size_t n, uint64_t* const h)
{
    for(std::size_t i = 0; i < n; i += L)
    {
        const auto c = std::min(L, n - i);  // Number of k-mers in the lanes.
        uint64_t x_l[L] = {}, h_l[L];
        std::copy(x + i, x + i + c, x_l);
        vec f;
        std::memcpy(&f, x_l, sizeof(f));

        // The reverse complements: complemented bases, reversed by swapping 2-bit pairs, then nibbles, then bytes.
        vec r = ~f;
        r = ((r >> 2) & splat(0x3333333333333333)) | ((r & splat(0x3333333333333333)) << 2);
        r = ((r >> 4) & splat(0x0F0F0F0F0F0F0F0F)) | ((r & splat(0x0F0F0F0F0F0F0F0F)) << 4);
        const auto byte_rev = _mm256_setr_epi8(7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8);
        r = reinterpret_cast<vec>(_mm256_shuffle_epi8(reinterpret_cast<__m256i>(r), byte_rev)) >> (64 - 2 * k);

        // The l-mer at offset `p` of a k-mer has its reverse complement at offset `k - l - p` of the k-mer's reverse
        // complement.
        vec m = splat(~uint64_t(0));
        for(std::size_t p = 0; p < w; ++p)
            m = min(m, wang_u64(min((f >> (2 * (k - l - p))) & splat(lmer_mask), (r >> (2 * p)) & splat(lmer_mask))));

        std::memcpy(h_l, &m, sizeof(m));
        std::copy(h_l, h_l + c, h + i);
    }
}


}



#endif
