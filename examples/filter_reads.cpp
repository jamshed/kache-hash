
// Filters reads by their share of k-mers in a reference, as K2Rmini does,
// though in either strand rather than in the reference's forward strand alone:
// builds the set of a reference's 31-mers from a number of threads, then
// prints the reads of a FASTA or FASTQ file that have at least the share
// `threshold` of their k-mer positions in the set, as FASTA, in their order in
// the input; and how many it kept, to stderr.
// Usage: filter_reads <reference.fa> <reads.fa|reads.fq> [threshold] [threads]

#include "fastx.hpp"

#include <kache_hash/kache_hash.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string>
#include <thread>
#include <vector>


constexpr uint16_t k = 31;


int main(const int argc, const char* const argv[])
{
    if(argc < 3)
    {
        std::cerr << "usage: filter_reads <reference.fa> <reads.fa|reads.fq> [threshold] [threads]\n";
        return 1;
    }

    const double threshold = argc > 3 ? std::stod(argv[3]) : 0.5;
    const std::size_t thread_c = argc > 4 ? std::stoul(argv[4]) : 4;
    if(!(threshold > 0 && threshold <= 1) || thread_c == 0)
    {
        std::cerr << "threshold: in (0, 1]; threads: at least 1\n";
        return 1;
    }

    const auto ref = examples::read_fastx(argv[1]);
    const auto reads = examples::read_fastx(argv[2]);

    // The reference's k-mers, inserted from a thread per part of its pieces.
    const auto pieces = examples::pieces(ref, k);
    kache_hash::Kmer_Set<k> set(std::max<std::size_t>(examples::positions(ref, k), 1));
    const std::size_t part = (pieces.size() + thread_c - 1) / thread_c;
    std::vector<std::thread> threads;
    for(std::size_t t = 0; t < thread_c; ++t)
        threads.emplace_back([&, t]()
        {
            const auto token = set.register_user();
            for(std::size_t i = t * part; i < std::min(pieces.size(), (t + 1) * part); ++i)
                set.insert(pieces[i].data(), pieces[i].size(), token);
        });
    for(auto& t : threads)
        t.join();

    // Whether each read stays: its k-mers in the set, counted from a thread per
    // part of the reads, as queries take no tokens.
    std::vector<char> keep(reads.size());
    const std::size_t read_part = (reads.size() + thread_c - 1) / thread_c;
    threads.clear();
    for(std::size_t t = 0; t < thread_c; ++t)
        threads.emplace_back([&, t]()
        {
            for(std::size_t r = t * read_part; r < std::min(reads.size(), (t + 1) * read_part); ++r)
            {
                const auto& s = reads[r].seq;
                const std::size_t positions = s.size() >= k ? s.size() - k + 1 : 0;
                keep[r] = positions > 0 && set.count(s.data(), s.size()) >= threshold * positions;
            }
        });
    for(auto& t : threads)
        t.join();

    std::size_t kept = 0;
    for(std::size_t r = 0; r < reads.size(); ++r)
        if(keep[r])
            std::cout << ">" << reads[r].name << "\n" << reads[r].seq << "\n", kept++;

    std::cerr << "kept " << kept << " of " << reads.size() << " reads, with at least " << threshold << " of their " << k << "-mers in the reference\n";
    return 0;
}
