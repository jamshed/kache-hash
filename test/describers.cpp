// Checks the k-mer describers, the batch and the SIMD lanes, against a brute-force reference: for every k-mer of the
// inputs, its canonical form, checksum, and minimizer hash, the least Wang hash over the canonical forms of all its
// l-mers. The inputs go whole, split at N, and cut into ragged pieces, so that stretches end at every offset of a
// chunk; the lanes take them 4 at a time. Both describers also go once with the k-mers' orientations, which the
// reference checks too. Also checks the minimizer hashes that resizes recompute from the canonical k-mers, in batches
// of every size up to a bucket's.
// Usage: describers <path>...

#include "kache_hash/Kmer_Batch.hpp"
#include "kache_hash/Kmer_Lanes.hpp"
#include "kache_hash/Kmer.hpp"
#include "Parser.hpp"
#include "kache_hash/DNA_Utility.hpp"

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>
#include <type_traits>


using namespace kache_hash;

constexpr uint16_t k = 31;
constexpr uint16_t l = 22;


// Returns the maximal stretches of at least `k` ACGT bases of the sequences at `path`.
static std::vector<std::string> read(const std::string& path)
{
    std::vector<std::string> F;
    Parser input(path);
    while(input.read_next_seq())
    {
        const char* const s = input.seq();
        const std::size_t len = input.seq_len();
        std::size_t start = 0;
        for(std::size_t i = 0; i <= len; ++i)
            if(i == len || DNA_Utility::is_placeholder(s[i]))
            {
                if(i - start >= k)
                    F.emplace_back(s + start, i - start);
                start = i + 1;
            }
    }

    input.close();
    return F;
}


// Returns the 2-bit encoding of the `n` bases at `s`, and of their reverse complement at `rc`.
static uint64_t encode(const char* const s, const std::size_t n, uint64_t& rc)
{
    uint64_t x = 0;
    rc = 0;
    for(std::size_t i = 0; i < n; ++i)
    {
        const uint64_t b = DNA_Utility::map_base(s[i]);
        x = (x << 2) | b;
        rc |= (3 - b) << (2 * i);
    }

    return x;
}


// Returns the minimizer hash of the k-mer at `s`: the least Wang hash of the canonical forms of its l-mers.
static uint64_t min_hash(const char* const s)
{
    uint64_t h = ~uint64_t(0);
    for(std::size_t o = 0; o + l <= k; ++o)
    {
        uint64_t rc;
        const auto x = encode(s + o, l, rc);
        h = std::min(h, wang_u64(std::min(x, rc)));
    }

    return h;
}


// Mismatch counts of the descriptions of some k-mers against the reference.
struct Tally
{
    std::size_t kmer_c = 0, key = 0, c = 0, h = 0, o = 0;

    // Checks the description `key`, `h` and `c` of the k-mer at `s`, and its orientation `fwd` unless it is negative:
    // whether the k-mer, as it occurs, is its canonical form.
    void check(const char* const s, const uint64_t key, const uint64_t h, const uint8_t c, const int fwd = -1)
    {
        uint64_t rc;
        const auto x = encode(s, k, rc);
        const auto canon = std::min(x, rc);
        kmer_c++;
        this->key += (key != canon || Kmer<k>(s).canonical().as_int() != canon);
        this->c += (c != kmer_checksum(canon) || c == 0);
        this->h += (h != min_hash(s));
        this->o += (fwd >= 0 && fwd != (x == canon));
    }

    bool bad() const { return key || c || h || o; }
};


// Checks both describers on the sequences `S` of the input `path`, which are
// `pieces` of it or its stretches, with the k-mers' orientations iff
// `orient_`. Returns the number of describers that failed.
template <bool orient_>
static int check(const char* const path, const std::vector<std::string>& S, const bool pieces)
{
    const char* const what = (pieces ? "pieces" : "stretches");
    const char* const how = (orient_ ? ", oriented" : "");

    Tally t;
    std::size_t expected = 0;
    Kmer_Batch<k, l> kb;
    std::conditional_t<orient_, typename Kmer_Batch<k, l>::Oriented_Chunk, typename Kmer_Batch<k, l>::Chunk> ch;
    for(const auto& s : S)
    {
        const auto n = s.size() - k + 1;
        expected += n;
        kb.reset(s.data(), n);
        for(std::size_t i = 0; kb.next(ch); i += ch.n)
            for(std::size_t j = 0; j < ch.n; ++j)
                if constexpr(orient_)
                    t.check(s.data() + i + j, ch.key[j], ch.h[j], ch.c[j], ch.o[j]);
                else
                    t.check(s.data() + i + j, ch.key[j], ch.h[j], ch.c[j]);
    }

    const bool bad = (t.bad() || t.kmer_c != expected);
    std::printf("describers %s, batch over %zu %s%s: k-mers %zu (ref %zu), bad keys %zu, checksums %zu, hashes %zu, orientations %zu%s\n",
                path, S.size(), what, how, t.kmer_c, expected, t.key, t.c, t.h, t.o, bad ? "  <-- WRONG" : "  ok");

    typedef Kmer_Lanes<k, l> lanes_t;
    constexpr auto L = lanes_t::L;
    Tally u;
    lanes_t kl;
    typename lanes_t::Chunk lc;
    for(std::size_t g = 0; g < S.size(); g += L)
    {
        const char* T[L];
        std::size_t len[L];
        const auto c = std::min(L, S.size() - g);
        for(std::size_t j = 0; j < c; ++j)
            T[j] = S[g + j].data(), len[j] = S[g + j].size();

        kl.reset(T, len, c);
        while(kl.template next<orient_>(lc))
            for(std::size_t j = 0; j < c; ++j)
            {
                const auto [t_b, t_e] = lc.range(j);
                for(std::size_t t = t_b; t < t_e; ++t)
                {
                    // The orientation, if asked for, is bit 8 of the checksum's lane; the bits above it stay clear.
                    const auto c_l = lc.c[t][j];
                    u.check(T[j] + lc.first + static_cast<int64_t>(t), lc.key[t][j], lc.h[t][j], static_cast<uint8_t>(c_l), orient_ ? int((c_l >> 8) & 1) : -1);
                    u.c += ((c_l >> (orient_ ? 9 : 8)) != 0);
                }
            }
    }

    const bool l_bad = (u.bad() || u.kmer_c != expected);
    std::printf("describers %s, lanes over %zu %s%s: k-mers %zu (ref %zu), bad keys %zu, checksums %zu, hashes %zu, orientations %zu%s\n",
                path, S.size(), what, how, u.kmer_c, expected, u.key, u.c, u.h, u.o, l_bad ? "  <-- WRONG" : "  ok");

    return bad + l_bad;
}


// Checks the minimizer hashes recomputed from the canonical k-mers of the
// sequences `S` of the input `path`, in batches whose sizes cycle from 1 to 32.
// Returns whether it failed.
static int check_recompute(const char* const path, const std::vector<std::string>& S)
{
    constexpr std::size_t max_batch = 32;
    uint64_t x[max_batch], r[max_batch], h[max_batch];  // A batch's canonical k-mers, and their reference and
                                                        // recomputed hashes.
    std::size_t c = 0, n = 1;   // Numbers of k-mers in the batch, and to fill it with.
    std::size_t kmer_c = 0, bad = 0;
    const auto flush = [&]()
    {
        Kmer_Lanes<k, l>::minimizer_hashes(x, c, h);
        for(std::size_t j = 0; j < c; ++j)
            bad += (h[j] != r[j]);

        kmer_c += c, c = 0, n = n % max_batch + 1;
    };

    for(const auto& s : S)
        for(std::size_t i = 0; i + k <= s.size(); ++i)
        {
            uint64_t rc;
            const auto f = encode(s.data() + i, k, rc);
            x[c] = std::min(f, rc), r[c] = min_hash(s.data() + i);
            if(++c == n)
                flush();
        }

    flush();

    std::printf("describers %s, minimizers recomputed over %zu stretches: k-mers %zu, bad hashes %zu%s\n",
                path, S.size(), kmer_c, bad, bad ? "  <-- WRONG" : "  ok");
    return bad > 0;
}


int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        std::fprintf(stderr, "Usage: %s <path>...\n", argv[0]);
        return 2;
    }

    int fail = 0;
    for(int a = 1; a < argc; ++a)
    {
        const auto F = read(argv[a]);
        std::vector<std::string> G; // Ragged pieces, overlapping by `k - 1`.
        for(const auto& s : F)
            for(std::size_t i = 0, len = 31 + (s.size() % 97); i + k <= s.size(); i += len - (k - 1), len = 31 + (len * 7 + 13) % 300)
                G.emplace_back(s.substr(i, std::min(len, s.size() - i)));

        const std::vector<std::string>* const inputs[] = {&F, &G};
        for(const auto* S : inputs)
            fail += check<false>(argv[a], *S, S == &G) + check<true>(argv[a], *S, S == &G);

        fail += check_recompute(argv[a], F);
    }

    return fail;
}
