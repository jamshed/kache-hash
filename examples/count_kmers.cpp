
// Counts the 31-mers of a FASTA or FASTQ file from a number of threads, and
// prints the number of distinct k-mers, of their occurrences, and of the
// k-mers by their counts, up to 10, as the Rust crate's example `kmer_counts`
// does.
// Usage: count_kmers <seqs.fa|seqs.fq> [threads]

#include "fastx.hpp"

#include <kache_hash/kache_hash.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string>
#include <thread>
#include <vector>


constexpr uint16_t k = 31;
constexpr std::size_t max_count = 10;   // The counts that the histogram tells apart: 1 to 9, and 10 or more.


int main(const int argc, const char* const argv[])
{
    if(argc < 2)
    {
        std::cerr << "usage: count_kmers <seqs.fa|seqs.fq> [threads]\n";
        return 1;
    }

    const std::size_t thread_c = argc > 2 ? std::stoul(argv[2]) : 4;
    if(thread_c == 0)
    {
        std::cerr << "threads: at least 1\n";
        return 1;
    }

    const auto recs = examples::read_fastx(argv[1]);
    const auto pieces = examples::pieces(recs, k);

    // The table grows beyond the k-mers it is sized for as needed; the k-mer
    // positions bound them, and so spare it the resizes.
    kache_hash::Kmer_Map<k, uint32_t> counts(std::max<std::size_t>(examples::positions(recs, k), 1));

    // A thread per part of the pieces, each registered for its own token.
    const std::size_t part = (pieces.size() + thread_c - 1) / thread_c;
    std::vector<std::thread> threads;
    for(std::size_t t = 0; t < thread_c; ++t)
        threads.emplace_back([&, t]()
        {
            const auto token = counts.register_user();
            for(std::size_t i = t * part; i < std::min(pieces.size(), (t + 1) * part); ++i)
                counts.upsert(pieces[i].data(), pieces[i].size(), [](const uint32_t c){ return c + 1; }, 1, token);
        });
    for(auto& t : threads)
        t.join();

    // The histogram and the occurrences, from a share of the table per thread.
    std::vector<std::array<uint64_t, max_count + 1>> hists(thread_c);
    std::vector<uint64_t> occurrences(thread_c);
    threads.clear();
    for(std::size_t t = 0; t < thread_c; ++t)
        threads.emplace_back([&, t]()
        {
            std::array<uint64_t, max_count + 1> h{};
            uint64_t o = 0;
            counts.for_each([&](const kache_hash::Kmer<k>&, const uint32_t c){ h[std::min<std::size_t>(c, max_count)]++, o += c; }, thread_c, t);
            hists[t] = h, occurrences[t] = o;
        });
    for(auto& t : threads)
        t.join();

    std::array<uint64_t, max_count + 1> hist{};
    uint64_t occurrence_c = 0;
    for(std::size_t t = 0; t < thread_c; ++t)
    {
        for(std::size_t c = 0; c <= max_count; ++c)
            hist[c] += hists[t][c];
        occurrence_c += occurrences[t];
    }

    std::cout << "distinct " << k << "-mers\t" << counts.size() << "\n";
    std::cout << "occurrences\t" << occurrence_c << "\n";
    for(std::size_t c = 1; c <= max_count; ++c)
        std::cout << "count " << c << (c == max_count ? "+" : "") << "\t" << hist[c] << "\n";

    return 0;
}
