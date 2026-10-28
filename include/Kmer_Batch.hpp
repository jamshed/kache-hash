
#ifndef KACHE_HASH_KMER_BATCH_HPP
#define KACHE_HASH_KMER_BATCH_HPP



#include "DNA.hpp"

#include <cstdint>
#include <cstddef>
#include <bit>
#include <algorithm>
#include <cassert>


namespace kache_hash
{


// Returns Thomas Wang's 64-bit integer hash of `x`: shifts, adds, and xors
// only, which SIMD lanes run without lowering the clock. It is a bijection on
// 64-bit values. The minimizers and the checksums come from it.
template <typename V_ = uint64_t>
inline constexpr V_ wang_u64(V_ x)
{
    x = (~x) + (x << 21);
    x = x ^ (x >> 24);
    x = (x + (x << 3)) + (x << 8);
    x = x ^ (x >> 14);
    x = (x + (x << 2)) + (x << 4);
    x = x ^ (x >> 28);
    x = x + (x << 31);
    return x;
}


// Returns the 8-bit checksum of the k-mer with the checksum hash `x`, never 0:
// its top byte. Its collision rate among the k-mers of a bucket matches a
// uniform byte's.
inline constexpr uint8_t checksum_of_hash(const uint64_t x)
{
    const auto c = static_cast<uint8_t>(x >> 56);
    return c ? c : 1;   // Checksum 0 marks empty slots.
}


// Returns the 8-bit checksum of the canonical k-mer with the 64-bit encoding
// `x`.
inline constexpr uint8_t kmer_checksum(const uint64_t x) { return checksum_of_hash(wang_u64(x)); }


// The minimizer of a k-mer is its l-mer with the least hash of the canonical
// form. Its coordinate in the k-mer is its offset in the strand its canonical
// form comes from, and that strand's orientation bit.
constexpr uint8_t coord_fwd_mask = 0b0100'0000;  // Bit of a minimizer coordinate set iff the minimizer is from the forward strand.
constexpr uint8_t coord_idx_mask = 0b0011'1111;  // Bits of a minimizer coordinate with the minimizer's offset.


// Describes the consecutive `k`-mers of a stretch of ACGT bases in chunks, for
// the streaming operations of the tables: their canonical forms, checksums,
// `l`-minimizer hashes, and minimizer coordinates. The minimizer of a k-mer is
// its l-mer of the least hash, so it is the least of the per-position minima;
// and the l-mer hash is a bijection, so any position holding the least hash has
// the minimizer. On ties, i.e. on a repeated minimizer, the leftmost position
// is taken.
template <uint16_t k, uint16_t l>
class Kmer_Batch
{
    static_assert(k <= 32, "Batched k-mer descriptions need k-mers in a machine word.");
    static_assert(l <= k && l <= 32);

public:

    static constexpr std::size_t C = 32;    // Maximum number of k-mers in a chunk.

private:

    static constexpr uint16_t w = k - l + 1;    // Number of l-mers in a k-mer.
    static constexpr std::size_t P = C + w - 1; // Maximum number of l-mer positions spanned by a chunk.
    static constexpr uint64_t kmer_mask = (k == 32 ? ~uint64_t(0) : (uint64_t(1) << (2 * k)) - 1);
    static constexpr uint64_t lmer_mask = (l == 32 ? ~uint64_t(0) : (uint64_t(1) << (2 * l)) - 1);

    static constexpr uint16_t lvl = std::bit_width(w) - 1;  // Number of doubling levels of the sliding minimum: `2^lvl` is the largest span within `w`.
    static constexpr std::size_t span = std::size_t(1) << lvl;  // That span.
    static constexpr std::size_t off = w - span;    // A window `[i, i + w)` is the union of the spans at `i` and `i + off`.

    const char* t;  // The stretch.
    std::size_t n;  // Number of k-mers in the stretch.
    std::size_t p;  // Index of the first k-mer of the next chunk.

    uint64_t fwd;   // Forward form of the last k-mer described, or of its prefix before the first.
    uint64_t rev;   // Reverse complement of the same.
    uint64_t l_fwd; // Forward form of the last l-mer hashed, or of its prefix before the first.
    uint64_t l_rev; // Reverse complement of the same.

    alignas(64) uint64_t hc[P]; // `hc[x]` is the least hash over the strands of the l-mer at the chunk's `x`'th position.
    alignas(64) uint64_t o_rev[P];  // `o_rev[x]` is 1 iff that least hash is of the reverse strand, else 0.
    std::size_t min_pos;    // Position of the last k-mer's minimizer, relative to the chunk, plus `w`; under `w` if outside the chunk.

    // Maps the ACGT character `ch`, in either case, to its 2-bit encoding.
    static constexpr uint64_t base(const char ch) { return ((ch >> 2) ^ (ch >> 1)) & 0b11; }

    // Puts the minima of the windows of `w` positions at `h[0..cnt)`, from the
    // minima `top` of the spans of `span` positions.
    static void window_min(const uint64_t* top, uint64_t* h, std::size_t cnt);

    // Puts the minima of the spans of `2 * s_` positions of `in` at `out`,
    // which holds the minima of the spans of `s_` positions over `pos_c`
    // positions, then recurses up to the span `span`.
    template <std::size_t s_> static void doubling(const uint64_t* in, uint64_t* out, std::size_t pos_c);

public:

    // Descriptions of a chunk of k-mers.
    struct Chunk
    {
        uint64_t key[C];    // `key[i]` is the canonical form of the chunk's `i`'th k-mer.
        uint64_t h[C];  // `h[i]` is its minimizer hash.
        uint8_t c[C];   // `c[i]` is its checksum, never 0.
        uint8_t m[C];   // `m[i]` is its minimizer coordinate, oriented per the canonical form, if described.
        std::size_t n;  // Number of k-mers in the chunk.
    };

    // Resets the describer to the stretch `t` of `n` k-mers.
    void reset(const char* t, std::size_t n);

    // Describes the next chunk of the stretch into `out`, and returns its
    // number of k-mers; 0 iff the stretch is exhausted. Describes minimizer
    // coordinates iff `coords_`, which is to be the same for all the chunks.
    template <bool coords_ = true> std::size_t next(Chunk& out);
};


template <uint16_t k, uint16_t l>
inline void Kmer_Batch<k, l>::reset(const char* const t, const std::size_t n)
{
    this->t = t, this->n = n, p = 0;

    fwd = rev = 0;
    for(std::size_t j = 0; j + 1 < k; ++j)
    {
        const auto b = base(t[j]);
        fwd = (fwd << 2) | b;
        rev = (rev >> 2) | ((DNA::T - b) << (2 * (k - 1)));
    }

    l_fwd = l_rev = 0;
    for(std::size_t j = 0; j + 1 < l; ++j)
    {
        const auto b = base(t[j]);
        l_fwd = (l_fwd << 2) | b;
        l_rev = (l_rev >> 2) | ((DNA::T - b) << (2 * (l - 1)));
    }

    min_pos = 0;
}


template <uint16_t k, uint16_t l>
template <std::size_t s_>
inline void Kmer_Batch<k, l>::doubling(const uint64_t* const in, uint64_t* const out, const std::size_t pos_c)
{
    if constexpr(s_ < span)
    {
        for(std::size_t x = 0; x + 2 * s_ <= pos_c; ++x)
            out[x] = std::min(in[x], in[x + s_]);

        doubling<2 * s_>(out, out + P, pos_c);
    }
}


template <uint16_t k, uint16_t l>
inline void Kmer_Batch<k, l>::window_min(const uint64_t* const top, uint64_t* const h, const std::size_t cnt)
{
    for(std::size_t i = 0; i < cnt; ++i)
        h[i] = std::min(top[i], top[i + off]);
}


template <uint16_t k, uint16_t l>
template <bool coords_>
inline std::size_t Kmer_Batch<k, l>::next(Chunk& out)
{
    out.n = 0;
    if(p >= n)
        return 0;

    const auto cnt = std::min(C, n - p);    // Number of k-mers in the chunk.
    const auto pos_c = cnt + w - 1;         // Number of l-mer positions they span; the first `w - 1` are kept from the last chunk, if any.

    // Hash the new l-mers.
    alignas(64) uint64_t lf[P], lr[P];
    const std::size_t x_0 = (p == 0 ? 0 : w - 1);
    for(std::size_t x = x_0; x < pos_c; ++x)
    {
        const auto b = base(t[p + x + l - 1]);
        l_fwd = ((l_fwd << 2) | b) & lmer_mask;
        l_rev = (l_rev >> 2) | ((DNA::T - b) << (2 * (l - 1)));
        lf[x] = l_fwd, lr[x] = l_rev;
    }

    for(std::size_t x = x_0; x < pos_c; ++x)
    {
        hc[x] = wang_u64(std::min(lf[x], lr[x]));
        if constexpr(coords_)
            o_rev[x] = (lr[x] < lf[x]);
    }

    // Sliding minima over `w` positions, by doubling spans; levels `j > 0` sit at `v + (j - 1) * P`.
    alignas(64) uint64_t v[lvl > 0 ? lvl * P : 1];
    const uint64_t* top = hc;   // Minima of the spans of `span` positions.
    if constexpr(lvl > 0)
    {
        doubling<1>(hc, v, pos_c);
        top = v + (lvl - 1) * P;
    }

    window_min(top, out.h, cnt);

    // The minimizers' positions. A new minimum can only enter at a window's right end; else the last k-mer's
    // minimizer stays, unless it has left the window, and only then is the window scanned, for its leftmost position
    // holding the least hash. The hash is a bijection, so the l-mer at such a position is the minimizer.
    if constexpr(coords_)
    {
        std::size_t pos = min_pos;  // Position of the last k-mer's minimizer, offset by `w` to keep it unsigned.
        for(std::size_t i = 0; i < cnt; ++i)
        {
            const auto h = out.h[i];
            if(hc[i + w - 1] == h)
                pos = i + w - 1 + w;
            else if(pos < i + w)    // The last minimizer has left the window.
            {
                uint8_t o = w - 1;
                for(std::size_t j = w - 1; j-- > 0; )
                    o = (hc[i + j] == h ? uint8_t(j) : o);
                pos = i + o + w;
            }

            const uint8_t o = pos - w - i;  // The minimizer's offset in the k-mer.
            const auto r = o_rev[pos - w];
            out.m[i] = uint8_t(r * (k - l - o)) | uint8_t((1 - r) * (o | coord_fwd_mask));
        }

        min_pos = (pos >= cnt + w ? pos - cnt : 0);   // Rebased to the next chunk, whose position 0 is this one's `cnt`.
    }

    // The k-mers themselves.
    for(std::size_t i = 0; i < cnt; ++i)
    {
        const auto b = base(t[p + i + k - 1]);
        fwd = ((fwd << 2) | b) & kmer_mask;
        rev = (rev >> 2) | ((DNA::T - b) << (2 * (k - 1)));

        const uint64_t fwd_mask = -uint64_t(fwd < rev); // All ones iff the forward form is canonical.
        out.key[i] = (fwd & fwd_mask) | (rev & ~fwd_mask);
        if constexpr(coords_)
            out.m[i] ^= (coord_fwd_mask & ~uint8_t(fwd_mask));
        out.c[i] = kmer_checksum(out.key[i]);
    }

    // Keep the last `w - 1` l-mer positions, the first of the next chunk.
    for(std::size_t x = 0; x + 1 < w; ++x)
    {
        hc[x] = hc[cnt + x];
        if constexpr(coords_)
            o_rev[x] = o_rev[cnt + x];
    }

    p += cnt;
    out.n = cnt;
    return cnt;
}


}



#endif
