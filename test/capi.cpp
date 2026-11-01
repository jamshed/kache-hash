
// Checks of the C interface, through it alone, against an exact reference. First argument validation and the bound on
// users; then a set, and maps of each width built by each operation, plain and oriented, presized and from a tiny table
// that resizes, and updated by each operation. The inserts and updates whose results don't depend on their order run
// from `threads` threads, the others from one. Each table is checked through iterations over its shares, its size,
// reductions and statistics, and per-sequence counts, masked counts and sums, per-position lookups and single lookups
// of queries: pieces of the index, their reverse complements, the negatives, and empty, short and all-N sequences.
// The inserts are of pieces of a part of the index, some reverse-complemented, some with lowercase runs and non-ACGT
// characters, in ragged batches. Last, forward-strand membership over the whole index, as K2Rmini uses it.
// Usage: capi <index> <negative> <threads>

#include "kache_capi.h"
#include "Parser.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>


using namespace kache_hash;

constexpr uint32_t k = 31;
constexpr uint64_t kmer_mask = (uint64_t(1) << (2 * k)) - 1;


// Returns the 2-bit encoding of the ACGT character `c`, in either case.
static uint64_t base(const char c) { return ((c >> 2) ^ (c >> 1)) & 3; }

static bool is_acgt(const char c) { const char x = c | 0x20; return x == 'a' || x == 'c' || x == 'g' || x == 't'; }


// Calls `f(i, key, fwd)` for each k-mer of `s` of ACGT characters only, at index `i`, with its canonical encoding `key`
// and whether it occurs as is.
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


// Calls `f(x)` with the forward encoding `x` of each k-mer of `s` of ACGT characters only.
template <typename F>
static void forward_kmers_of(const std::string& s, const F& f)
{
    uint64_t fwd = 0;
    std::size_t run = 0;
    for(const char c : s)
        if(!is_acgt(c))
            run = 0;
        else if(fwd = ((fwd << 2) | base(c)) & kmer_mask; ++run >= k)
            f(fwd);
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


// Returns pieces of the sequences at `path` of at most 4096 bases, overlapping by `k - 1`.
static std::vector<std::string> pieces(const std::string& path)
{
    std::vector<std::string> S;
    Parser input(path);
    while(input.read_next_seq())
    {
        const std::string s(input.seq(), input.seq_len());
        for(std::size_t i = 0; i + k <= s.size(); i += 4096 - (k - 1))
            S.emplace_back(s.substr(i, 4096));
    }

    input.close();
    return S;
}


// Returns the number of k-mer positions of the sequences `S[b..e)`.
static uint64_t positions(const std::vector<std::string>& S, const std::size_t b, const std::size_t e)
{
    uint64_t p = 0;
    for(std::size_t i = b; i < e; ++i)
        p += (S[i].size() >= k ? S[i].size() - k + 1 : 0);

    return p;
}


// Returns `op(v, x)`, for values of at most `max`.
static uint64_t apply(const uint32_t op, const uint64_t v, const uint64_t x, const uint64_t max)
{
    switch(op)
    {
    case KACHE_KEEP:    return v;
    case KACHE_ASSIGN:  return x;
    case KACHE_ADD:     return v > max - x ? max : v + x;
    case KACHE_OR:      return v | x;
    case KACHE_MIN:     return std::min(v, x);
    default:            return std::max(v, x);
    }
}


// Returns the operand of the operation `op` for the call `b`, for values of `bits` bits: some sums saturate, and the
// ORs set one bit each.
static uint64_t operand(const uint32_t op, const uint64_t b, const uint32_t bits)
{
    const uint64_t max = (bits == 64 ? ~uint64_t(0) : (uint64_t(1) << bits) - 1);
    const uint64_t h = (b + 1) * 0x9E3779B97F4A7C15;
    switch(op)
    {
    case KACHE_ADD: return b % 5 == 0 ? max / 2 + 1 : 1 + b % 7;
    case KACHE_OR:  return uint64_t(1) << ((h >> 32) % bits);
    default:        return (h >> 7) & max;
    }
}


// Returns whether the results of the operation `op` don't depend on the order of its applications.
static bool commutes(const uint32_t op) { return op != KACHE_KEEP && op != KACHE_ASSIGN; }


// Calls `f(b, i, c)` for each of the batches `b` of the sequences `S`, of `c` sequences from `S[i]` each, from `T`
// threads, each taking every `T`'th batch; or from this thread, in order, if `T` is 0. Batches are of 1 to 13
// sequences.
template <typename F>
static void batches(const std::vector<std::string>& S, const std::size_t T, const F& f)
{
    std::vector<std::pair<std::size_t, std::size_t>> B;  // The batches' first sequences and sizes.
    for(std::size_t i = 0, b = 0; i < S.size(); i += B.back().second, ++b)
        B.emplace_back(i, std::min<std::size_t>(1 + b % 13, S.size() - i));

    if(T == 0)
    {
        for(std::size_t b = 0; b < B.size(); ++b)
            f(0, b, B[b].first, B[b].second);
        return;
    }

    std::vector<std::thread> W;
    for(std::size_t t = 0; t < T; ++t)
        W.emplace_back([&, t]()
        {
            for(std::size_t b = t; b < B.size(); b += T)
                f(t, b, B[b].first, B[b].second);
        });
    for(auto& w : W)
        w.join();
}


// Puts the pointers and lengths of the sequences `S[i..i + c)` at `seqs` and `lens`.
static void spans(const std::vector<std::string>& S, const std::size_t i, const std::size_t c, std::vector<const char*>& seqs, std::vector<size_t>& lens)
{
    seqs.resize(c), lens.resize(c);
    for(std::size_t j = 0; j < c; ++j)
        seqs[j] = S[i + j].data(), lens[j] = S[i + j].size();
}


// Keys and values of a table, from an iteration.
struct Iteration
{
    uint32_t width;
    std::vector<std::pair<uint64_t, uint64_t>> kv;
    std::size_t chunks = 0;
    bool stop = false;  // Whether to ask for no more chunks after the first.
};

static int collect(void* const ctx, const uint64_t* const keys, const void* const values, const size_t n)
{
    auto& it = *static_cast<Iteration*>(ctx);
    for(size_t i = 0; i < n; ++i)
    {
        uint64_t v = 0;
        if(values != nullptr)
            std::memcpy(&v, static_cast<const uint8_t*>(values) + i * it.width, it.width);
        it.kv.emplace_back(keys[i], v);
    }

    it.chunks++;
    return it.stop;
}


// A table's expected contents: each key's value, and for keys whose last oriented assignment had them in both
// orientations, the other value they may have.
struct Reference
{
    std::unordered_map<uint64_t, uint64_t> val;
    std::unordered_map<uint64_t, uint64_t> alt;
};


// Checks the table `t` of values of `width` bytes against its reference `R`, through iterations over `T` shares from
// `T` threads, its size, reductions and statistics; adopts the table's values of the keys with alternatives. Returns
// the number of keys and values wrong.
static std::size_t check_contents(const kache_table* const t, const uint32_t width, Reference& R, const std::size_t T, std::size_t& errors)
{
    std::vector<Iteration> I(T, Iteration{width, {}});
    std::vector<int> st(T);
    std::vector<std::thread> W;
    for(std::size_t s = 0; s < T; ++s)
        W.emplace_back([&, s](){ st[s] = kache_for_each(t, T, s, collect, &I[s]); });
    for(auto& w : W)
        w.join();

    std::vector<std::pair<uint64_t, uint64_t>> kv;
    for(std::size_t s = 0; s < T; ++s)
    {
        errors += (st[s] != KACHE_OK);
        kv.insert(kv.end(), I[s].kv.begin(), I[s].kv.end());
    }

    std::sort(kv.begin(), kv.end());
    std::size_t wrong = (kv.size() != R.val.size());
    for(std::size_t i = 0; i < kv.size(); ++i)
    {
        const auto [x, v] = kv[i];
        const auto it = R.val.find(x);
        if((i > 0 && kv[i - 1].first == x) || it == R.val.end())
        {
            wrong++;
            continue;
        }

        if(v != it->second)
        {
            const auto a = R.alt.find(x);
            if(a != R.alt.end() && a->second == v)
                it->second = v;
            else
                wrong++;
        }
    }

    // A stopped iteration receives one chunk.
    Iteration one{width, {}, 0, true};
    errors += (kache_for_each(t, 1, 0, collect, &one) != KACHE_OK || one.chunks != (R.val.empty() ? 0 : 1));

    uint64_t sz = 0;
    errors += (kache_size(t, &sz) != KACHE_OK || sz != R.val.size());

    const uint64_t max = (width == 8 ? ~uint64_t(0) : (uint64_t(1) << (8 * width)) - 1);
    for(const uint64_t m : {max, uint64_t(0x5A5A5A5A5A5A5A5A) & max})
        for(uint32_t what = KACHE_REDUCE_KEYS; what <= uint32_t(width == 0 ? KACHE_REDUCE_KEYS : KACHE_REDUCE_SUM); ++what)
        {
            uint64_t ref = 0;
            for(const auto& [x, v] : R.val)
                ref += (width == 0 ? 1 : what == KACHE_REDUCE_KEYS ? ((v & m) != 0) : what == KACHE_REDUCE_POPCOUNT ? std::popcount(v & m) : (v & m));

            uint64_t r = 0;
            errors += (kache_reduce(t, what, m, T, &r) != KACHE_OK || r != ref);
        }

    kache_stats stats;
    errors += (kache_get_stats(t, &stats) != KACHE_OK || stats.keys != R.val.size() || stats.capacity < stats.keys ||
               (stats.keys > 0 && std::abs(stats.primary_pct + stats.secondary_pct + stats.overflow_pct - 100) > 1e-6));

    return wrong;
}


// Expected results of a query sequence.
struct Expected
{
    uint64_t count = 0, masked = 0, sum = 0;
    std::vector<uint8_t> found, fwd;
    std::vector<uint64_t> val;
};

// Returns the expected results of the query `s` in a table with the reference `R`, for the masks `m_fwd` and `m_rev`.
static Expected expected(const std::string& s, const Reference& R, const uint64_t m_fwd, const uint64_t m_rev)
{
    Expected e;
    const auto p = (s.size() >= k ? s.size() - k + 1 : 0);
    e.found.assign(p, KACHE_NOT_KMER), e.fwd.assign(p, 0), e.val.assign(p, 0);
    kmers_of(s, [&](const std::size_t i, const uint64_t x, const bool f)
    {
        const auto it = R.val.find(x);
        e.found[i] = (it != R.val.end() ? KACHE_PRESENT : KACHE_ABSENT), e.fwd[i] = f;
        if(it != R.val.end())
        {
            e.count++, e.sum += it->second, e.val[i] = it->second;
            e.masked += ((it->second & (f ? m_fwd : m_rev)) != 0);
        }
    });

    return e;
}


// Checks the queries `Q` of the table `t` of values of `width` bytes, with the reference `R`, from `T` threads: their
// counts, masked counts, sums, lookups, and single lookups. Returns the number of sequences with wrong results.
static std::size_t check_queries(const kache_table* const t, const uint32_t width, const Reference& R, const std::vector<std::string>& Q, const std::size_t T, std::size_t& errors)
{
    const uint64_t m_fwd = 0x5555555555555555, m_rev = 0xA3A3A3A3A3A3A3A3;
    std::vector<std::size_t> wrong(T, 0), err(T, 0);
    batches(Q, T, [&](const std::size_t th, const std::size_t b, const std::size_t i, const std::size_t c)
    {
        std::vector<const char*> seqs;
        std::vector<size_t> lens;
        spans(Q, i, c, seqs, lens);
        std::vector<uint64_t> cnt(c), msk(c), sum(c);
        const auto p = positions(Q, i, i + c);
        std::vector<uint8_t> found(p), fwd(p), val(p * width);
        const bool want_val = (width > 0 && b % 2 == 0), want_fwd = (b % 3 != 1);   // Each lookup variant, some calls each.

        err[th] += (kache_count(t, seqs.data(), lens.data(), c, cnt.data()) != KACHE_OK);
        err[th] += (kache_lookup(t, seqs.data(), lens.data(), c, found.data(), want_val ? val.data() : nullptr, want_fwd ? fwd.data() : nullptr) != KACHE_OK);
        if(width > 0)
        {
            err[th] += (kache_count_masked(t, seqs.data(), lens.data(), c, m_fwd, m_rev, msk.data()) != KACHE_OK);
            err[th] += (kache_sum(t, seqs.data(), lens.data(), c, sum.data()) != KACHE_OK);
        }

        for(std::size_t j = 0, o = 0; j < c; ++j)
        {
            const auto e = expected(Q[i + j], R, m_fwd, m_rev);
            bool bad = (cnt[j] != e.count) || (width > 0 && (msk[j] != e.masked || sum[j] != e.sum));
            for(std::size_t x = 0; x < e.found.size(); ++x, ++o)
            {
                uint64_t v = 0;
                if(want_val)
                    std::memcpy(&v, val.data() + o * width, width);
                bad |= (found[o] != e.found[x]) || (want_fwd && fwd[o] != e.fwd[x]) || (want_val && v != e.val[x]);
            }

            // The sequence's first k-mer, alone.
            if(Q[i + j].size() >= k)
            {
                uint8_t f = 2;
                uint64_t v = 1;
                err[th] += (kache_get(t, Q[i + j].data(), &f, &v) != KACHE_OK);
                const bool has = (e.found[0] == KACHE_PRESENT);
                bad |= (f != has) || (v != (has ? e.val[0] : 0));
            }

            wrong[th] += bad;
        }
    });

    std::size_t w = 0;
    for(std::size_t th = 0; th < T; ++th)
        w += wrong[th], errors += err[th];

    return w;
}


int main(const int argc, const char* const argv[])
{
    if(argc < 4)
    {
        std::fprintf(stderr, "Usage: %s <index> <negative> <threads>\n", argv[0]);
        return 2;
    }

    const std::size_t T = std::clamp<std::size_t>(std::stoul(argv[3]), 1, KACHE_MAX_USERS);
    int fail = 0;
    const auto expect = [&](const bool ok, const char* const what)
    {
        if(!ok)
            std::printf("capi: %s  FAILED\n", what), fail++;
    };

    // Argument validation.
    kache_table* t = nullptr;
    expect(kache_supports(31, 0) && kache_supports(31, 1) && kache_supports(31, 2) && kache_supports(31, 4) && kache_supports(31, 8), "supports");
    expect(!kache_supports(25, 0) && !kache_supports(31, 3) && !kache_supports(31, 16), "doesn't support");
    auto o = kache_default_options(31, 4);
    expect(o.k == 31 && o.value_bytes == 4 && o.expected_keys == (1 << 22) && o.load_factor == 0.8 && o.resize_workers == 16, "default options");
    expect(kache_new(&o, nullptr) == KACHE_EINVAL && kache_new(nullptr, &t) == KACHE_EINVAL && t == nullptr, "new rejects null arguments");
    for(const auto& [k_, w_] : {std::pair<uint32_t, uint32_t>{25, 4}, {31, 3}})
    {
        auto x = o;
        x.k = k_, x.value_bytes = w_;
        expect(kache_new(&x, &t) == KACHE_EINVAL && t == nullptr, "new rejects unsupported tables");
    }
    for(const double lf : {0.0, 1.0, 1.5, std::nan("")})
    {
        auto x = o;
        x.load_factor = lf;
        expect(kache_new(&x, &t) == KACHE_EINVAL && t == nullptr, "new rejects the load factor");
    }
    o.resize_workers = 0;
    expect(kache_new(&o, &t) == KACHE_EINVAL && t == nullptr, "new rejects 0 workers");
    kache_free(nullptr);

    // The bound on users, and arguments of the operations.
    o = kache_default_options(31, 1), o.expected_keys = 1000, o.resize_workers = 2;
    expect(kache_new(&o, &t) == KACHE_OK && t != nullptr, "new");
    uint32_t u = 0;
    bool reg_ok = true;
    for(uint32_t i = 0; i < KACHE_MAX_USERS; ++i)
        reg_ok &= (kache_register(t, &u) == KACHE_OK && u == i);
    expect(reg_ok, "128 registrations");
    expect(kache_register(t, &u) == KACHE_EUSERS, "the 129th registration is refused");

    const char* const one = "ATGCGTACCTTAGGCATCCAGTTAGACGGATCA";  // 3 distinct k-mers.
    const size_t one_len = 33, zero_len = 0, five = 5;
    const char* const null_seq = nullptr;
    uint64_t c = 0;
    uint8_t f = 0, b[8];
    expect(kache_insert(t, KACHE_MAX_USERS, &one, &one_len, 1, KACHE_ADD, 1, nullptr) == KACHE_EINVAL, "an unregistered user is refused");
    expect(kache_insert(t, 0, &one, &one_len, 1, KACHE_MAX + 1, 1, nullptr) == KACHE_EINVAL, "an unknown operation is refused");
    expect(kache_insert(t, 0, &one, &one_len, 1, KACHE_ADD, 256, nullptr) == KACHE_EINVAL, "an operand too wide is refused");
    expect(kache_insert(t, 0, nullptr, &one_len, 1, KACHE_ADD, 1, nullptr) == KACHE_EINVAL, "null sequences are refused");
    expect(kache_insert(t, 0, &null_seq, &five, 1, KACHE_ADD, 1, nullptr) == KACHE_EINVAL, "a null sequence is refused");
    expect(kache_insert(t, 0, &null_seq, &zero_len, 1, KACHE_ADD, 1, &c) == KACHE_OK && c == 0, "a null empty sequence is taken");
    for(const uint32_t op : {uint32_t(KACHE_KEEP), uint32_t(KACHE_MIN)})
        expect(kache_insert_oriented(t, 0, &one, &one_len, 1, op, 1, 2, nullptr) == KACHE_EINVAL, "an oriented insert refuses the operation");
    expect(kache_insert_oriented(t, 0, &one, &one_len, 1, KACHE_OR, 1, 256, nullptr) == KACHE_EINVAL, "an oriented insert refuses an operand too wide");
    expect(kache_insert(t, 0, &one, &one_len, 1, KACHE_ADD, 255, &c) == KACHE_OK && c == 3, "insert");
    expect(kache_update(t, &one, &one_len, 1, KACHE_ADD, 1, &c) == KACHE_OK && c == 3, "update");
    expect(kache_get(t, one, &f, &c) == KACHE_OK && f == 1 && c == 255, "updates saturate");
    expect(kache_get(t, one, nullptr, &c) == KACHE_EINVAL, "get refuses a null output");
    expect(kache_count(t, &one, &one_len, 1, nullptr) == KACHE_EINVAL, "count refuses a null output");
    expect(kache_lookup(t, &one, &one_len, 1, nullptr, nullptr, nullptr) == KACHE_EINVAL, "lookup refuses a null output");
    expect(kache_reduce(t, KACHE_REDUCE_SUM + 1, ~uint64_t(0), 1, &c) == KACHE_EINVAL && kache_reduce(t, KACHE_REDUCE_SUM, 1, 0, &c) == KACHE_EINVAL, "reduce refuses its arguments");
    expect(kache_for_each(t, 2, 2, collect, nullptr) == KACHE_EINVAL && kache_for_each(t, 1, 0, nullptr, nullptr) == KACHE_EINVAL, "for_each refuses its arguments");
    kache_free(t), t = nullptr;

    o = kache_default_options(31, 0), o.expected_keys = 1000, o.resize_workers = 2;
    expect(kache_new(&o, &t) == KACHE_OK && kache_register(t, &u) == KACHE_OK, "new set");
    expect(kache_insert(t, u, &one, &one_len, 1, KACHE_MAX + 7, ~uint64_t(0), &c) == KACHE_OK && c == 3, "a set ignores the operation");
    expect(kache_insert_oriented(t, u, &one, &one_len, 1, KACHE_OR, 1, 2, nullptr) == KACHE_EINVAL, "a set refuses oriented inserts");
    expect(kache_update(t, &one, &one_len, 1, KACHE_ADD, 1, nullptr) == KACHE_EINVAL, "a set refuses updates");
    expect(kache_count_masked(t, &one, &one_len, 1, 1, 2, &c) == KACHE_EINVAL && kache_sum(t, &one, &one_len, 1, &c) == KACHE_EINVAL, "a set refuses value queries");
    expect(kache_lookup(t, &one, &one_len, 1, b, b, nullptr) == KACHE_EINVAL, "a set's lookups have no values");
    expect(kache_reduce(t, KACHE_REDUCE_POPCOUNT, 1, 1, &c) == KACHE_EINVAL && kache_reduce(t, KACHE_REDUCE_KEYS, 0, 1, &c) == KACHE_OK && c == 3, "a set reduces its keys");
    kache_free(t), t = nullptr;

    // The inserts: a part of the index, with every third piece reverse-complemented, the reverse complement of every
    // seventh beside it, and lowercase runs and non-ACGT characters in every eleventh.
    const auto index = pieces(argv[1]);
    std::vector<std::string> S;
    for(std::size_t i = 0; i < index.size() && positions(S, 0, S.size()) < 1'000'000; ++i)
    {
        auto s = (i % 3 == 0 ? revcomp(index[i]) : index[i]);
        if(i % 11 == 5 && s.size() > 300)
        {
            std::transform(s.begin() + 100, s.begin() + 160, s.begin() + 100, [](const char x){ return x | 0x20; });
            s[40] = 'N', s[200] = 'n', s[230] = 'R', s[260] = '-';
        }

        S.push_back(s);
        if(i % 7 == 0)
            S.push_back(revcomp(s));
    }

    // The queries: every fourth piece, the reverse complement of every ninth, the negatives, and some edge cases.
    std::vector<std::string> Q;
    for(std::size_t i = 0; i < S.size(); ++i)
    {
        if(i % 4 == 0)
            Q.push_back(S[i]);
        if(i % 9 == 0)
            Q.push_back(revcomp(S[i]));
    }
    const auto neg = pieces(argv[2]);
    Q.insert(Q.end(), neg.begin(), neg.begin() + std::min<std::size_t>(neg.size(), 150));
    Q.insert(Q.end(), {"", "ACGT", std::string(40, 'N'), S[0].substr(0, k), S[0].substr(0, k - 1), "acgtnACGTACGTACGTACGTACGTACGTACGTACGTAC"});

    struct Scenario { uint32_t width, op; bool oriented; uint32_t update_op; };
    std::vector<Scenario> X{{0, KACHE_KEEP, false, KACHE_MAX + 1}, {0, KACHE_KEEP, false, KACHE_MAX + 1}};
    for(const uint32_t w : {1, 2, 4, 8})
    {
        for(const uint32_t op : {KACHE_KEEP, KACHE_ASSIGN, KACHE_ADD, KACHE_OR, KACHE_MIN, KACHE_MAX})
            X.push_back({w, op, false, op == KACHE_KEEP ? uint32_t(KACHE_ASSIGN) : op == KACHE_ASSIGN ? uint32_t(KACHE_KEEP) : op});
        for(const uint32_t op : {KACHE_ASSIGN, KACHE_ADD, KACHE_OR, KACHE_MAX})
            X.push_back({w, op, true, KACHE_MAX + 1});
    }

    const char* const op_name[] = {"keep", "assign", "add", "or", "min", "max"};
    for(std::size_t s = 0; s < X.size(); ++s)
    {
        const auto& x = X[s];
        const uint32_t bits = 8 * x.width;
        const uint64_t max = (x.width == 8 ? ~uint64_t(0) : (uint64_t(1) << bits) - 1);
        const bool tiny = (s % 2 == 1);
        const std::size_t threads = (commutes(x.op) || x.width == 0 ? T : 0);
        std::size_t errors = 0;

        // The reference, by the inserts in their order.
        Reference R;
        std::unordered_map<uint64_t, std::pair<std::size_t, uint8_t>> last;  // An oriented assignment's last call of each key, and its orientations there.
        batches(S, 0, [&](std::size_t, const std::size_t b, const std::size_t i, const std::size_t c)
        {
            const auto x_fwd = operand(x.op, b, std::max(bits, 1u)), x_rev = operand(x.op, b + 977, std::max(bits, 1u));
            for(std::size_t j = i; j < i + c; ++j)
                kmers_of(S[j], [&](std::size_t, const uint64_t y, const bool fwd)
                {
                    const auto v = (x.width == 0 ? 0 : !x.oriented ? x_fwd : fwd ? x_fwd : x_rev);
                    const auto it = R.val.find(y);
                    if(it == R.val.end())
                        R.val.emplace(y, x.oriented ? apply(x.op, 0, v, max) : v);
                    else if(x.width > 0)
                        it->second = apply(x.op, it->second, v, max);

                    if(x.oriented && x.op == KACHE_ASSIGN)
                    {
                        auto& [l_b, l_o] = last[y];
                        if(l_b != b || l_o == 0)
                            l_b = b, l_o = 0;
                        l_o |= (fwd ? 1 : 2);
                    }
                });
        });
        for(const auto& [y, lb] : last)
            if(lb.second == 3)  // Either orientation's operand may come last.
            {
                const auto x_fwd = operand(x.op, lb.first, bits), x_rev = operand(x.op, lb.first + 977, bits);
                R.alt[y] = (R.val[y] == x_fwd ? x_rev : x_fwd);
            }

        // The table, presized or tiny, by the same inserts.
        auto opt = kache_default_options(31, x.width);
        opt.expected_keys = (tiny ? 10 : R.val.size()), opt.resize_workers = std::max<std::size_t>(T, 1);
        expect(kache_new(&opt, &t) == KACHE_OK, "new");
        std::vector<uint32_t> user(std::max<std::size_t>(threads, 1));
        for(auto& y : user)
            errors += (kache_register(t, &y) != KACHE_OK);

        std::vector<uint64_t> ins(user.size(), 0);
        std::vector<std::size_t> err(user.size(), 0);
        batches(S, threads, [&](const std::size_t th, const std::size_t b, const std::size_t i, const std::size_t c)
        {
            std::vector<const char*> seqs;
            std::vector<size_t> lens;
            spans(S, i, c, seqs, lens);
            uint64_t n = 0;
            const auto x_fwd = operand(x.op, b, std::max(bits, 1u)), x_rev = operand(x.op, b + 977, std::max(bits, 1u));
            const auto r = (x.oriented ? kache_insert_oriented(t, user[th], seqs.data(), lens.data(), c, x.op, x_fwd, x_rev, &n)
                                       : kache_insert(t, user[th], seqs.data(), lens.data(), c, x.op, x.width > 0 ? x_fwd : 0, &n));
            err[th] += (r != KACHE_OK), ins[th] += n;
        });
        uint64_t inserted = 0;
        for(std::size_t th = 0; th < user.size(); ++th)
            inserted += ins[th], errors += err[th];

        const auto wrong_kv = check_contents(t, x.width, R, T, errors);
        const auto wrong_q = check_queries(t, x.width, R, Q, T, errors);

        // The updates, by the queries, then the contents again.
        std::size_t wrong_up = 0, wrong_present = 0;
        if(x.update_op <= KACHE_MAX)
        {
            std::vector<std::size_t> bad(T, 0), e(T, 0);
            const auto up_threads = (commutes(x.update_op) ? T : 0);
            batches(Q, up_threads, [&](const std::size_t th, const std::size_t b, const std::size_t i, const std::size_t c)
            {
                std::vector<const char*> seqs;
                std::vector<size_t> lens;
                spans(Q, i, c, seqs, lens);
                uint64_t present = 0, ref = 0;
                e[th] += (kache_update(t, seqs.data(), lens.data(), c, x.update_op, operand(x.update_op, b + 5000, bits), &present) != KACHE_OK);
                for(std::size_t j = i; j < i + c; ++j)
                    kmers_of(Q[j], [&](std::size_t, const uint64_t y, bool){ ref += R.val.count(y); });
                bad[th] += (present != ref);
            });
            batches(Q, 0, [&](std::size_t, const std::size_t b, const std::size_t i, const std::size_t c)
            {
                const auto v = operand(x.update_op, b + 5000, bits);
                for(std::size_t j = i; j < i + c; ++j)
                    kmers_of(Q[j], [&](std::size_t, const uint64_t y, bool){ if(const auto it = R.val.find(y); it != R.val.end()) it->second = apply(x.update_op, it->second, v, max); });
            });
            for(std::size_t th = 0; th < T; ++th)
                wrong_present += bad[th], errors += e[th];

            wrong_up = check_contents(t, x.width, R, T, errors);
        }

        kache_stats st;
        kache_get_stats(t, &st);
        kache_free(t), t = nullptr;

        const bool ok = (errors == 0 && inserted == R.val.size() && wrong_kv == 0 && wrong_q == 0 && wrong_up == 0 && wrong_present == 0);
        fail += !ok;
        char name[64];
        std::snprintf(name, sizeof(name), "%s%s %s", x.width == 0 ? "set" : "map", x.width == 0 ? "" : std::to_string(bits).c_str(), x.width == 0 ? "insert" : op_name[x.op]);
        std::printf("capi %s %s%s%s %s %s: keys %zu, ins %llu, placed %.2f/%.2f/%.2f%%, contents %zu wrong, queries %zu of %zu wrong",
                    argv[1], name, x.oriented ? " oriented" : "", x.update_op <= KACHE_MAX ? (std::string(", then ") + op_name[x.update_op]).c_str() : "",
                    threads > 0 ? ("T=" + std::to_string(threads)).c_str() : "in order", tiny ? "tiny" : "presized", R.val.size(), (unsigned long long)inserted,
                    st.primary_pct, st.secondary_pct, st.overflow_pct, wrong_kv, wrong_q, Q.size());
        if(x.update_op <= KACHE_MAX)
            std::printf(", updates %zu presence and %zu contents wrong", wrong_present, wrong_up);
        std::printf(", %zu errors  %s\n", errors, ok ? "ok" : "FAILED");
    }

    // Forward-strand membership, as K2Rmini uses it, over the whole index from a tiny table: every third piece goes in
    // reverse-complemented; a k-mer is a member iff inserted as it is.
    std::vector<std::string> P = index;
    for(std::size_t i = 0; i < P.size(); i += 3)
        P[i] = revcomp(P[i]);

    std::vector<uint64_t> F;    // The distinct forward k-mers inserted.
    for(const auto& s : P)
        forward_kmers_of(s, [&](const uint64_t y){ F.push_back(y); });
    std::sort(F.begin(), F.end());
    F.erase(std::unique(F.begin(), F.end()), F.end());

    auto opt = kache_default_options(31, 1);
    opt.expected_keys = 10, opt.resize_workers = std::max<std::size_t>(T, 1);
    expect(kache_new(&opt, &t) == KACHE_OK, "new, tiny");
    std::vector<uint32_t> user(T);
    for(auto& y : user)
        expect(kache_register(t, &y) == KACHE_OK, "register");
    std::vector<int> st(T, KACHE_OK);
    batches(P, T, [&](const std::size_t th, std::size_t, const std::size_t i, const std::size_t c)
    {
        std::vector<const char*> seqs;
        std::vector<size_t> lens;
        spans(P, i, c, seqs, lens);
        if(const auto r = kache_insert_oriented(t, user[th], seqs.data(), lens.data(), c, KACHE_OR, 1, 2, nullptr); r != KACHE_OK)
            st[th] = r;
    });
    expect(std::all_of(st.begin(), st.end(), [](const int x){ return x == KACHE_OK; }), "inserts");

    uint64_t members = 0;
    expect(kache_reduce(t, KACHE_REDUCE_POPCOUNT, 3, T, &members) == KACHE_OK && members == F.size(), "members");

    // Per-sequence counts of members, from `T` threads: every piece, the reverse complement of every fifth, and the
    // negatives.
    std::vector<std::string> M = P;
    for(std::size_t i = 0; i < P.size(); i += 5)
        M.push_back(revcomp(P[i]));
    M.insert(M.end(), neg.begin(), neg.end());
    M.insert(M.end(), {"", "ACGT", std::string(40, 'N')});
    std::vector<std::size_t> wrong(T, 0);
    batches(M, T, [&](const std::size_t th, std::size_t, const std::size_t i, const std::size_t c)
    {
        std::vector<const char*> seqs;
        std::vector<size_t> lens;
        spans(M, i, c, seqs, lens);
        std::vector<uint64_t> counts(c);
        if(kache_count_masked(t, seqs.data(), lens.data(), c, 1, 2, counts.data()) != KACHE_OK)
            return (void)(wrong[th] += c);

        for(std::size_t j = 0; j < c; ++j)
        {
            uint64_t ref = 0;
            forward_kmers_of(M[i + j], [&](const uint64_t y){ ref += std::binary_search(F.begin(), F.end(), y); });
            wrong[th] += (counts[j] != ref);
        }
    });
    std::size_t wrong_c = 0;
    for(const auto w : wrong)
        wrong_c += w;
    expect(wrong_c == 0, "per-sequence member counts");

    kache_stats stats;
    expect(kache_get_stats(t, &stats) == KACHE_OK && stats.keys > 0 && stats.keys <= F.size(), "stats");
    std::printf("capi %s forward-strand T=%zu: members %llu (ref %zu), keys %llu, capacity %llu, placed %.2f/%.2f/%.2f%%, counts %zu wrong over %zu sequences  %s\n",
                argv[1], T, (unsigned long long)members, F.size(), (unsigned long long)stats.keys, (unsigned long long)stats.capacity,
                stats.primary_pct, stats.secondary_pct, stats.overflow_pct, wrong_c, M.size(), fail ? "FAILED" : "ok");
    kache_free(t);
    return fail;
}
