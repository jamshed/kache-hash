// Regression test for concurrent insertion into `Streaming_Kmer_Hash_Table`.
//
// Several threads insert the same k-mers, in the same order, into one table
// (a hash set, then a hash map; created for the distinct k-mer count, then
// grown from a small capacity so that concurrent resizes happen too). Running
// in lockstep, the threads contend on each bucket at the moment its first
// k-mer is placed. After each round the test checks that the table's size
// equals the number of distinct k-mers and that every k-mer is found.
//
// Usage: concurrent_insert_test [sequence-length [rounds [threads]]]
//        (defaults: 4'000'000 bases, 20 rounds, min(4, hardware threads)).
// Exit status: 0 if every round is correct, 1 otherwise.

#include "Streaming_Kmer_Hash_Table.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <thread>
#include <vector>


constexpr uint16_t k = 31;  // k-mer size.
constexpr uint16_t l = 19;  // Minimizer size (the table's default).


// Calls `f(w)` on each k-mer window of the sequence `s`.
template <typename F>
void for_each_kmer(const std::string& s, const F& f)
{
    kache_hash::Kmer_Window<k, l> w;
    w.init(s.data());
    for(std::size_t i = 0; ; ++i)
    {
        f(w);
        if(i + k == s.size())
            break;

        w.advance(s[i + k]);
    }
}


// Inserts the k-mers of `seqs` into a table of type `Table` from `thread_c`
// threads, `rounds` times, and returns the number of incorrect rounds.
// `max_sz` is the initial capacity of the table: the distinct k-mer count, or
// smaller to trigger resizes.
template <typename Table, bool is_map>
std::size_t run(const std::string& name, const std::vector<std::string>& seqs, const std::size_t distinct, const std::size_t max_sz, const std::size_t rounds, const std::size_t thread_c)
{
    std::size_t bad = 0;
    for(std::size_t r = 0; r < rounds; ++r)
    {
        Table T(max_sz, thread_c);
        std::vector<typename Table::Token> tokens;
        for(std::size_t t = 0; t < thread_c; ++t)
            tokens.push_back(T.register_user());

        std::vector<std::thread> workers;
        for(std::size_t t = 0; t < thread_c; ++t)
            workers.emplace_back([&, t]()
            {
                for(const auto& s : seqs)
                    for_each_kmer(s, [&](const auto& w)
                    {
                        if constexpr(is_map)
                            T.insert(w, uint64_t(1), tokens[t]);
                        else
                            T.insert(w, tokens[t]);
                    });
            });

        for(auto& w : workers)
            w.join();

        std::size_t missing = 0;
        for(const auto& s : seqs)
            for_each_kmer(s, [&](const auto& w){ missing += !T.find(w); });

        if(T.size() != distinct || missing > 0)
        {
            bad++;
            std::cerr << name << ", round " << r << ": size " << T.size() << ", expected " << distinct
                      << "; " << missing << " k-mer lookups failed.\n";
        }
    }

    std::cout << name << ": " << (rounds - bad) << " of " << rounds << " rounds correct.\n";
    return bad;
}


int main(int argc, char* argv[])
{
    const std::size_t len = argc > 1 ? std::size_t(std::strtoull(argv[1], nullptr, 10)) : 4'000'000;
    const std::size_t rounds = argc > 2 ? std::size_t(std::strtoull(argv[2], nullptr, 10)) : 20;
    const std::size_t hw = std::max(1u, std::thread::hardware_concurrency());
    // The table keeps per-hardware-thread state indexed by the token id.
    const std::size_t thread_c = std::min<std::size_t>(hw, argc > 3 ? std::strtoull(argv[3], nullptr, 10) : std::min<std::size_t>(4, hw));

    // A random sequence, cut into overlapping pieces of 1000 bases.
    std::mt19937_64 rng(42);
    std::string g(len, 'A');
    for(auto& c : g)
        c = "ACGT"[rng() & 3];

    std::vector<std::string> seqs;
    for(std::size_t p = 0; p + k <= len; p += 1000)
        seqs.push_back(g.substr(p, std::min<std::size_t>(1000 + k - 1, len - p)));

    // Distinct k-mer count, from the single-threaded table.
    kache_hash::Streaming_Kmer_Hash_Table<k, false, void, l> S(len);
    const auto token = S.register_user();
    for(const auto& s : seqs)
        for_each_kmer(s, [&](const auto& w){ S.insert(w, token); });
    const std::size_t distinct = S.size();

    std::cout << len << " bases, " << distinct << " distinct k-mers, " << thread_c << " threads, " << rounds << " rounds.\n";

    using Set = kache_hash::Streaming_Kmer_Hash_Table<k, true, void, l>;
    using Map = kache_hash::Streaming_Kmer_Hash_Table<k, true, uint64_t, l>;
    std::size_t bad = 0;
    bad += run<Set, false>("hash set, sized for the k-mers", seqs, distinct, distinct, rounds, thread_c);
    bad += run<Map, true>("hash map, sized for the k-mers", seqs, distinct, distinct, rounds, thread_c);
    bad += run<Set, false>("hash set, growing (concurrent resizes)", seqs, distinct, 1 << 16, std::max<std::size_t>(rounds / 4, 1), thread_c);
    bad += run<Map, true>("hash map, growing (concurrent resizes)", seqs, distinct, 1 << 16, std::max<std::size_t>(rounds / 4, 1), thread_c);

    std::cout << (bad ? "FAILED" : "PASSED") << "\n";
    return bad ? 1 : 0;
}
