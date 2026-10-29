
// Correctness checks of the table against an exact reference: sizes, presence, counts, and negative hits, for sets and
// maps, sequential and concurrent tables, presized and resizing.
// Usage: check <index> <negative> <seq-set|seq-map|mt-set|mt-map>[8|16|32] <threads> <init-cap> <stream|raw|multi>
//   Maps hold 64-bit counts, or 8-, 16- or 32-bit ones per the suffix.
//   stream: operations over the fragments split at N.
//   raw:    operations over the raw sequences, N's included, which the table must split itself.
//   multi:  multi-sequence operations over batches of 13 raw pieces, of ragged lengths.

#include "Streaming_Kmer_Hash_Table.hpp"
#include "Parser.hpp"
#include "utility.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <thread>
#include <vector>


using namespace kache_hash;

constexpr uint16_t k = 31;
constexpr uint16_t l = 19;


// Returns the sequences at `path`: raw iff `raw`, else their maximal stretches of at least `k` ACGT bases.
static std::vector<std::string> read(const std::string& path, const bool raw)
{
    std::vector<std::string> F;
    Parser input(path);
    while(input.read_next_seq())
    {
        const char* const s = input.seq();
        const std::size_t len = input.seq_len();
        if(raw)
        {
            F.emplace_back(s, len);
            continue;
        }

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


static uint64_t canonical(const char* const s) { return Kmer<k>(s).canonical().as_int(); }


// Sorted distinct canonical k-mers of the fragments `F`, with their multiplicities.
struct Reference
{
    std::vector<uint64_t> kmer, count;

    explicit Reference(const std::vector<std::string>& F)
    {
        std::vector<uint64_t> K;
        for(const auto& s : F)
            for(std::size_t i = 0; i + k <= s.size(); ++i)
                K.push_back(canonical(s.data() + i));
        std::sort(K.begin(), K.end());
        for(std::size_t i = 0; i < K.size(); )
        {
            std::size_t j = i;
            while(j < K.size() && K[j] == K[i])
                j++;
            kmer.push_back(K[i]), count.push_back(j - i);
            i = j;
        }
    }

    std::size_t size() const { return kmer.size(); }
    uint64_t count_of(const uint64_t x) const { return count[std::lower_bound(kmer.begin(), kmer.end(), x) - kmer.begin()]; }
    bool has(const uint64_t x) const { return std::binary_search(kmer.begin(), kmer.end(), x); }
};


template <bool mt_, typename T_>
static int check(const std::string& index, const std::string& neg, const std::size_t T, const std::size_t init_cap, const std::string& mode)
{
    typedef Streaming_Kmer_Hash_Table<k, mt_, T_, l> ht_t;
    constexpr bool is_set = std::is_void_v<T_>;
    const auto inc = [](const uint64_t v){ return v + 1; };

    const auto F = read(index, false);
    std::vector<std::string> S; // What the table is fed: pieces of at most 4096 bases (ragged for `multi`), overlapping by `k - 1`, so that threads share sequences.
    for(const auto& s : (mode == "raw" || mode == "multi" ? read(index, true) : F))
        for(std::size_t i = 0, len = 4096; i + k <= s.size(); i += len - (k - 1), len = (mode == "multi" ? 31 + (len * 7 + 13) % 9000 : 4096))
            S.emplace_back(s.substr(i, len));
    const std::size_t batch_sz = (mode == "multi" ? 13 : 1);
    const Reference R(F);

    ht_t M(init_cap > 0 ? init_cap : R.size(), T);
    const auto cap_0 = M.capacity();
    std::vector<typename ht_t::Token> token(T);
    for(auto& t : token)
        t = M.register_user();

    std::vector<std::size_t> ins(T, 0);
    std::vector<std::thread> W;
    for(std::size_t t = 0; t < T; ++t)
        W.emplace_back([&, t]()
        {
            for(std::size_t i = t * batch_sz; i < S.size(); i += T * batch_sz)   // Pieces round-robin, so threads contend on the same buckets.
            {
                const char* seqs[13];
                std::size_t lens[13], c = 0;
                for(std::size_t j = i; j < std::min(S.size(), i + batch_sz); ++j)
                    seqs[c] = S[j].data(), lens[c] = S[j].size(), c++;

                if constexpr(is_set)
                    ins[t] += (mode == "multi" ? M.insert(seqs, lens, c, token[t]) : M.insert(seqs[0], lens[0], token[t]));
                else
                    ins[t] += (mode == "multi" ? M.upsert(seqs, lens, c, inc, 1, token[t]) : M.upsert(seqs[0], lens[0], inc, 1, token[t]));
            }
        });
    for(auto& w : W)
        w.join();

    std::size_t ins_c = 0;
    for(const auto x : ins)
        ins_c += x;

    // Lookups of every index k-mer, then of every negative k-mer; and counts of the negatives over the raw sequences.
    // All split over the threads.
    const auto NF = read(neg, false), NR = read(neg, true);
    const std::size_t L = std::max<std::size_t>(T, std::thread::hardware_concurrency());
    std::vector<std::size_t> c_miss(L), c_wrong(L), c_calls(L), c_neg(L), c_neg_wrong(L), c_neg_exp(L), c_neg_count(L);
    std::vector<std::thread> Q;
    for(std::size_t t = 0; t < L; ++t)
        Q.emplace_back([&, t]()
        {
            for(std::size_t i = t; i < F.size(); i += L)
            {
                const auto& s = F[i];
                const auto found = [&](const std::size_t q, const std::size_t j, const auto& v)
                {
                    c_calls[t]++;
                    if(q != 0 || !v)
                        c_miss[t]++;
                    else if constexpr(!is_set)  // Counts wrap at narrow widths.
                        c_wrong[t] += (*v != static_cast<T_>(R.count_of(canonical(s.data() + j))));
                };

                if(mode == "multi")
                {
                    const char* seqs[2] = {s.data(), s.data()};
                    const std::size_t lens[2] = {s.size(), 0};
                    M.find(seqs, lens, 2, found);
                }
                else
                    M.find(s.data(), s.size(), [&](const std::size_t j, const auto& v){ found(0, j, v); });
            }

            for(std::size_t i = t; i < NF.size(); i += L)
            {
                const auto& s = NF[i];
                M.find(s.data(), s.size(), [&](const std::size_t j, const auto& v)
                {
                    const bool exp = R.has(canonical(s.data() + j));
                    c_neg[t] += static_cast<bool>(v), c_neg_exp[t] += exp;
                    c_neg_wrong[t] += (static_cast<bool>(v) != exp);
                });
            }

            if(mode == "multi")
            {
                for(std::size_t i = t * 13; i < NR.size(); i += L * 13)
                {
                    const char* seqs[13];
                    std::size_t lens[13], c = 0;
                    for(std::size_t j = i; j < std::min(NR.size(), i + 13); ++j)
                        seqs[c] = NR[j].data(), lens[c] = NR[j].size(), c++;
                    c_neg_count[t] += M.count(seqs, lens, c);
                }
            }
            else
                for(std::size_t i = t; i < NR.size(); i += L)
                    c_neg_count[t] += M.count(NR[i].data(), NR[i].size());
        });
    for(auto& q : Q)
        q.join();

    const auto sum = [](const std::vector<std::size_t>& v){ std::size_t x = 0; for(const auto y : v) x += y; return x; };
    const auto miss = sum(c_miss), wrong = sum(c_wrong), calls = sum(c_calls);
    const auto neg_hits = sum(c_neg), neg_wrong = sum(c_neg_wrong), neg_expected = sum(c_neg_exp), neg_count = sum(c_neg_count);

    std::size_t kmer_c = 0;
    for(const auto& s : F)
        kmer_c += s.size() - k + 1;

    const bool bad = (M.size() != R.size() || ins_c != R.size() || miss || wrong || calls != kmer_c || neg_wrong || neg_count != neg_expected);
    std::printf("check %s %s-%s T=%zu cap %zu->%zu %s: size %zu ins %zu (ref %zu), miss %zu, wrong %zu, calls %zu (ref %zu), neg %zu/%zu (ref %zu, %zu wrong)%s\n",
                index.c_str(), mt_ ? "mt" : "seq", is_set ? "set" : "map", T, static_cast<std::size_t>(cap_0), static_cast<std::size_t>(M.capacity()), mode.c_str(),
                M.size(), ins_c, R.size(), miss, wrong, calls, kmer_c, neg_hits, neg_count, neg_expected, neg_wrong, bad ? "  <-- WRONG" : "  ok");
    return bad ? 1 : 0;
}


int main(int argc, char* argv[])
{
    if(argc != 7)
    {
        std::fprintf(stderr, "Usage: %s <index> <negative> <seq-set|seq-map|mt-set|mt-map>[8|16|32] <threads> <init-cap> <stream|raw|multi>\n", argv[0]);
        return 2;
    }

    const std::string index(argv[1]), neg(argv[2]), kind(argv[3]), mode(argv[6]);
    const std::size_t T = std::strtoull(argv[4], nullptr, 10), cap = std::strtoull(argv[5], nullptr, 10);
    if(mode != "stream" && mode != "raw" && mode != "multi")
        return 2;
    if(kind == "seq-set")   return check<false, void>(index, neg, 1, cap, mode);
    if(kind == "seq-map")   return check<false, uint64_t>(index, neg, 1, cap, mode);
    if(kind == "seq-map32") return check<false, uint32_t>(index, neg, 1, cap, mode);
    if(kind == "seq-map16") return check<false, uint16_t>(index, neg, 1, cap, mode);
    if(kind == "seq-map8")  return check<false, uint8_t>(index, neg, 1, cap, mode);
    if(kind == "mt-set")    return check<true, void>(index, neg, T, cap, mode);
    if(kind == "mt-map")    return check<true, uint64_t>(index, neg, T, cap, mode);
    if(kind == "mt-map32")  return check<true, uint32_t>(index, neg, T, cap, mode);
    if(kind == "mt-map16")  return check<true, uint16_t>(index, neg, T, cap, mode);
    if(kind == "mt-map8")   return check<true, uint8_t>(index, neg, T, cap, mode);
    return 2;
}
