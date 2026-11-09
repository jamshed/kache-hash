
// Correctness checks of the oriented operations against an exact reference: `upsert_oriented` with a mask of the
// orientations each k-mer occurs in, and with counts per orientation, then the oriented `find`'s orientations, hits and
// per-sequence counts, for concurrent tables, presized and resizing. The inputs are pieces of the index in either
// orientation, some with lowercase runs and non-ACGT characters, and groups of k-mers sharing a minimizer, which fill
// their buckets into the overflow table; the queries add the negatives, and empty, short and all-N sequences.
// Usage: oriented <index> <negative> <threads>

#include "kache_hash/Streaming_Kmer_Hash_Table.hpp"
#include "Parser.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <random>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>


using namespace kache_hash;

constexpr uint16_t k = 31;
constexpr uint16_t l = 22;

constexpr uint64_t kmer_mask = (uint64_t(1) << (2 * k)) - 1;


// Returns the 2-bit encoding of the ACGT character `c`, in either case, as the describers encode it.
static uint64_t base(const char c) { return ((c >> 2) ^ (c >> 1)) & 3; }

static bool is_acgt(const char c) { const char x = c | 0x20; return x == 'a' || x == 'c' || x == 'g' || x == 't'; }


// Calls `f(i, key, fwd)` for each k-mer of `s` of ACGT characters only, at index `i`, with its canonical encoding `key`
// and whether it occurs in that form.
template <typename F>
static void kmers_of(const std::string& s, const F& f)
{
    uint64_t fwd = 0, rev = 0;
    std::size_t run = 0;    // Number of ACGT characters ending at the current one.
    for(std::size_t j = 0; j < s.size(); ++j)
    {
        if(!is_acgt(s[j]))
        {
            run = 0;
            continue;
        }

        const auto b = base(s[j]);
        fwd = ((fwd << 2) | b) & kmer_mask;
        rev = (rev >> 2) | ((3 - b) << (2 * (k - 1)));
        if(++run >= k)
            f(j + 1 - k, std::min(fwd, rev), fwd <= rev);
    }
}


// Returns the reverse complement of `s`, keeping the case, and non-ACGT characters as they are.
static std::string revcomp(const std::string& s)
{
    std::string r(s.rbegin(), s.rend());
    for(auto& c : r)
        switch(c)
        {
        case 'A': c = 'T'; break;   case 'C': c = 'G'; break;   case 'G': c = 'C'; break;   case 'T': c = 'A'; break;
        case 'a': c = 't'; break;   case 'c': c = 'g'; break;   case 'g': c = 'c'; break;   case 't': c = 'a'; break;
        default: break;
        }

    return r;
}


// Returns the sequences at `path`, as they are.
static std::vector<std::string> read(const std::string& path)
{
    std::vector<std::string> F;
    Parser input(path);
    while(input.read_next_seq())
        F.emplace_back(input.seq(), input.seq_len());

    input.close();
    return F;
}


// Returns pieces of the sequences `F` of at most `len` bases, overlapping by `k - 1`.
static std::vector<std::string> pieces(const std::vector<std::string>& F, const std::size_t len)
{
    std::vector<std::string> S;
    for(const auto& s : F)
        for(std::size_t i = 0; i + k <= s.size(); i += len - (k - 1))
            S.emplace_back(s.substr(i, len));

    return S;
}


// Returns the minimizer hash of the l-mer at `s`, as the describers compute it.
static uint64_t lmer_hash(const char* const s)
{
    uint64_t f = 0, r = 0;
    for(std::size_t i = 0; i < l; ++i)
        f = (f << 2) | base(s[i]), r = r | ((3 - base(s[i])) << (2 * i));

    return wang_u64(std::min(f, r));
}

// Returns the minimizer hash of the k-mer at `s`: the least of its l-mers'.
static uint64_t kmer_min_hash(const char* const s)
{
    uint64_t h = UINT64_MAX;
    for(std::size_t i = 0; i + l <= k; ++i)
        h = std::min(h, lmer_hash(s + i));

    return h;
}

// Returns `G` groups of `n` distinct k-mers each, the k-mers of a group sharing a minimizer, and thus their primary and
// secondary buckets: only 96 of a group fit in the main table, if the groups' buckets don't collide.
static std::vector<std::vector<std::string>> minimizer_groups(const std::size_t G, const std::size_t n, const uint64_t seed)
{
    std::mt19937_64 rng(seed);
    const auto random_seq = [&](const std::size_t len){ std::string s(len, 'A'); for(auto& c : s) c = "ACGT"[rng() & 3]; return s; };

    std::vector<std::vector<std::string>> K(G);
    std::unordered_set<uint64_t> seen;
    for(auto& group : K)
    {
        std::string X;
        do  // An l-mer with a hash among the least 1/1024'th: likely the minimizer of k-mers containing it.
            X = random_seq(l);
        while(lmer_hash(X.data()) > (UINT64_MAX >> 10));

        const auto target = lmer_hash(X.data());
        while(group.size() < n)
        {
            auto s = random_seq(k);
            std::copy(X.begin(), X.end(), s.begin() + rng() % (k - l + 1));
            if(kmer_min_hash(s.data()) != target)
                continue;
            uint64_t key = 0;
            kmers_of(s, [&](std::size_t, const uint64_t x, bool){ key = x; });
            if(seen.insert(key).second)
                group.push_back(s);
        }
    }

    return K;
}


// The distinct canonical k-mers of what is inserted, with their occurrences in either orientation.
struct Reference
{
    std::vector<uint64_t> key;
    std::vector<uint32_t> fwd, rev;

    explicit Reference(const std::vector<std::string>& S)
    {
        std::vector<uint64_t> O;    // Occurrences: the key, shifted, with its orientation in the low bit.
        for(const auto& s : S)
            kmers_of(s, [&](std::size_t, const uint64_t x, const bool f){ O.push_back((x << 1) | f); });
        std::sort(O.begin(), O.end());
        for(std::size_t i = 0; i < O.size(); )
        {
            const auto x = O[i] >> 1;
            uint32_t f = 0, r = 0;
            for(; i < O.size() && (O[i] >> 1) == x; ++i)
                (O[i] & 1 ? f : r)++;
            key.push_back(x), fwd.push_back(f), rev.push_back(r);
        }
    }

    std::size_t size() const { return key.size(); }

    // Returns the index of the key `x`, or `size()` if absent.
    std::size_t at(const uint64_t x) const
    {
        const auto it = std::lower_bound(key.begin(), key.end(), x);
        return (it != key.end() && *it == x ? it - key.begin() : size());
    }

    uint8_t mask(const std::size_t i) const { return (fwd[i] > 0 ? 1 : 0) | (rev[i] > 0 ? 2 : 0); }

    // Returns whether the k-mer with the canonical encoding `x`, occurring in orientation `f`, was inserted so.
    bool has(const uint64_t x, const bool f) const { const auto i = at(x); return i < size() && (f ? fwd[i] : rev[i]) > 0; }
};


// Upserts the sequences `S` into the map `M` with their orientations, through `f`, from `T` threads in ragged batches,
// round-robin so that threads contend on the same buckets. Returns the number of k-mers inserted.
template <typename M_, typename F>
static std::size_t build(M_& M, const std::vector<std::string>& S, const std::size_t T, const F& f)
{
    std::vector<typename M_::Token> token(T);
    for(auto& t : token)
        t = M.register_user();

    std::vector<std::size_t> ins(T, 0);
    std::vector<std::thread> W;
    for(std::size_t t = 0; t < T; ++t)
        W.emplace_back([&, t]()
        {
            const char* seqs[13];
            std::size_t lens[13];
            for(std::size_t i = t * 13, b = 0; i < S.size(); i += T * 13, ++b) // A block of 13, a call each, or one at a time.
            {
                const std::size_t n = std::min<std::size_t>(13, S.size() - i), step = (b % 2 ? 13 : 1);
                for(std::size_t s = 0; s < n; s += step)
                {
                    const std::size_t c = std::min(step, n - s);
                    for(std::size_t j = 0; j < c; ++j)
                        seqs[j] = S[i + s + j].data(), lens[j] = S[i + s + j].size();
                    ins[t] += M.upsert_oriented(seqs, lens, c, f, token[t]);
                }
            }
        });

    for(auto& w : W)
        w.join();

    std::size_t n = 0;
    for(const auto x : ins)
        n += x;

    return n;
}


int main(const int argc, const char* const argv[])
{
    if(argc < 4)
    {
        std::fprintf(stderr, "Usage: %s <index> <negative> <threads>\n", argv[0]);
        return 2;
    }

    const std::string index = argv[1], neg = argv[2];
    const std::size_t T = std::stoul(argv[3]);

    // What is inserted: pieces of the index; every third reverse-complemented; every fourth with a lowercase run, and
    // every seventh with non-ACGT characters; then minimizer groups, inserted forward, reverse-complemented, or both.
    auto S = pieces(read(index), 4096);
    for(std::size_t i = 0; i < S.size(); ++i)
    {
        auto& s = S[i];
        if(i % 4 == 1)
            for(std::size_t j = s.size() / 4; j < s.size() / 2; ++j)
                s[j] |= 0x20;
        if(i % 7 == 2 && s.size() > 200)
            s[50] = 'N', s[90] = 'n', s[130] = 'R', s[170] = '-';
        if(i % 3 == 0)
            s = revcomp(s);
    }

    const auto groups = minimizer_groups(6, 160, 42);
    for(std::size_t g = 0; g < groups.size(); ++g)
        for(const auto& x : groups[g])
        {
            if(g % 3 != 1)
                S.push_back(x);
            if(g % 3 != 0)
                S.push_back(revcomp(x));
        }

    const Reference R(S);
    std::size_t fwd_c = 0;  // Number of distinct forward k-mers.
    for(std::size_t i = 0; i < R.size(); ++i)
        fwd_c += (R.fwd[i] > 0) + (R.rev[i] > 0);

    // The queries: what is inserted, its reverse complement for every fifth sequence, the negatives' pieces, and edge
    // cases, in ragged batches.
    std::vector<std::string> Q = S;
    for(std::size_t i = 0; i < S.size(); i += 5)
        Q.push_back(revcomp(S[i]));
    for(auto& s : pieces(read(neg), 4096))
        Q.push_back(std::move(s));
    Q.insert(Q.end(), {"", "ACGT", std::string(k - 1, 'A'), std::string(40, 'N'), std::string(k, 'a'), "NNNNACGTACGTACGTACGTACGTACGTACGTACGTNNNN"});

    int fail = 0;
    for(const std::size_t init_cap : {std::size_t(0), std::size_t(10)})
    {
        // The orientation mask.
        Streaming_Kmer_Hash_Table<k, true, uint8_t, l> M(init_cap > 0 ? init_cap : R.size(), T);
        const auto cap_0 = M.capacity();
        const auto ins = build(M, S, T, [](const uint8_t v, const bool fwd){ return static_cast<uint8_t>(v | (fwd ? 1 : 2)); });
        const auto [prim, sec, ovf] = M.key_distribution();

        std::size_t it_c = 0, it_wrong = 0, fwd_seen = 0;
        M.for_each([&](const Kmer<k>& x, const uint8_t v)
        {
            const auto i = R.at(x.as_int());
            it_c++, it_wrong += (i == R.size() || v != R.mask(i)), fwd_seen += std::popcount(v);
        }, 1, 0);

        // The oriented lookups, from `T` threads in ragged batches.
        std::vector<std::vector<uint8_t>> calls(Q.size());  // Callbacks per k-mer index of each query.
        for(std::size_t q = 0; q < Q.size(); ++q)
            calls[q].assign(Q[q].size() >= k ? Q[q].size() - k + 1 : 0, 0);
        std::vector<uint32_t> cnt(Q.size(), 0);
        std::vector<std::size_t> wrong_fwd(T, 0), wrong_hit(T, 0), wrong_cnt(T, 0), presence_wrong(T, 0);
        std::vector<std::thread> W;
        for(std::size_t t = 0; t < T; ++t)
            W.emplace_back([&, t]()
            {
                const char* seqs[13];
                std::size_t lens[13];
                uint32_t counts[13];
                for(std::size_t i_0 = t * 13, b = 0; i_0 < Q.size(); i_0 += T * 13, ++b)    // A block of 13, a call each, or 3 at a time.
                for(std::size_t i = i_0, n = std::min<std::size_t>(13, Q.size() - i_0), step = (b % 2 ? 13 : 3); i < i_0 + n; i += step)
                {
                    const std::size_t c = std::min(step, i_0 + n - i);
                    for(std::size_t j = 0; j < c; ++j)
                        seqs[j] = Q[i + j].data(), lens[j] = Q[i + j].size(), counts[j] = 0;

                    M.find(seqs, lens, c, [&](const std::size_t q, const std::size_t x, const std::optional<uint8_t>& r, const bool fwd)
                    {
                        uint64_t key = 0; bool f = false;
                        kmers_of(Q[i + q].substr(x, k), [&](std::size_t, const uint64_t y, const bool o){ key = y, f = o; });
                        calls[i + q][x]++;
                        wrong_fwd[t] += (fwd != f);
                        const bool hit = r && (*r & (fwd ? 1 : 2));
                        wrong_hit[t] += (hit != R.has(key, f));
                        counts[q] += hit;
                    });

                    M.find(seqs, lens, c, [&](const std::size_t q, const std::size_t x, const std::optional<uint8_t>& r)
                    {
                        uint64_t key = 0;
                        kmers_of(Q[i + q].substr(x, k), [&](std::size_t, const uint64_t y, bool){ key = y; });
                        presence_wrong[t] += (r.has_value() != (R.at(key) < R.size()));
                    });

                    for(std::size_t j = 0; j < c; ++j)
                    {
                        uint32_t ref = 0;
                        kmers_of(Q[i + j], [&](std::size_t, const uint64_t y, const bool o){ ref += R.has(y, o); });
                        wrong_cnt[t] += (counts[j] != ref);
                        cnt[i + j] = counts[j];
                    }
                }
            });
        for(auto& w : W)
            w.join();

        std::size_t call_wrong = 0, wf = 0, wh = 0, wc = 0, pw = 0;
        for(std::size_t q = 0; q < Q.size(); ++q)
        {
            std::vector<uint8_t> expect(calls[q].size(), 0);
            kmers_of(Q[q], [&](const std::size_t x, uint64_t, bool){ expect[x] = 1; });
            for(std::size_t x = 0; x < expect.size(); ++x)
                call_wrong += (calls[q][x] != expect[x]);
        }
        for(std::size_t t = 0; t < T; ++t)
            wf += wrong_fwd[t], wh += wrong_hit[t], wc += wrong_cnt[t], pw += presence_wrong[t];

        // The counts per orientation, which catch skipped, repeated or misoriented updates.
        Streaming_Kmer_Hash_Table<k, true, uint64_t, l> C(init_cap > 0 ? init_cap : R.size(), T);
        const auto ins_c = build(C, S, T, [](const uint64_t v, const bool fwd){ return v + (fwd ? 1 : uint64_t(1) << 32); });
        std::size_t cnt_c = 0, cnt_wrong = 0;
        C.for_each([&](const Kmer<k>& x, const uint64_t v)
        {
            const auto i = R.at(x.as_int());
            cnt_c++, cnt_wrong += (i == R.size() || (v & 0xFFFFFFFF) != R.fwd[i] || (v >> 32) != R.rev[i]);
        }, 1, 0);

        const bool spilled = (init_cap > 0 || (sec > 0 && M.overflow_size() > 0));  // The groups reach the overflow table when presized.
        const bool ok = ins == R.size() && M.size() == R.size() && it_c == R.size() && it_wrong == 0 && fwd_seen == fwd_c &&
                        call_wrong == 0 && wf == 0 && wh == 0 && wc == 0 && pw == 0 &&
                        ins_c == R.size() && C.size() == R.size() && cnt_c == R.size() && cnt_wrong == 0 && spilled;
        std::printf("oriented %s T=%zu cap %zu->%zu: keys %zu ins %zu (ref %zu), forward %zu (ref %zu), iterate %zu wrong, "
                    "placed %.2f/%.2f/%.2f%%, calls %zu wrong, fwd %zu wrong, hits %zu wrong, counts %zu wrong, presence %zu wrong, "
                    "per-orientation counts %zu/%zu wrong  %s\n",
                    index.c_str(), T, cap_0, M.capacity(), M.size(), ins, R.size(), fwd_seen, fwd_c, it_wrong,
                    prim, sec, ovf, call_wrong, wf, wh, wc, pw, cnt_wrong, cnt_c, ok ? "ok" : "FAILED");
        fail += !ok;
    }

    return fail;
}
