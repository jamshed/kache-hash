
// Tells, for each read of a FASTA or FASTQ file, how many of its 31-mers occur
// in a reference's forward strand, and how many in its reverse strand, through
// a map whose values flag the orientations each k-mer occurs in in the
// reference, as the Rust crate's example `strand_membership` does. Prints a
// line per read: its index, its number of k-mer positions, and its k-mers in
// each strand.
// Usage: strand_membership <reference.fa> <reads.fa|reads.fq> [threads]

#include "fastx.hpp"

#include <kache_hash/kache_hash.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <optional>
#include <string>
#include <thread>
#include <vector>


constexpr uint16_t k = 31;

// A k-mer's flags: it occurs in the reference as is, i.e. in its canonical
// form, or as its reverse complement.
constexpr uint8_t as_is = 1;
constexpr uint8_t revcomp = 2;


int main(const int argc, const char* const argv[])
{
    if(argc < 3)
    {
        std::cerr << "usage: strand_membership <reference.fa> <reads.fa|reads.fq> [threads]\n";
        return 1;
    }

    const std::size_t thread_c = argc > 3 ? std::stoul(argv[3]) : 4;
    if(thread_c == 0)
    {
        std::cerr << "threads: at least 1\n";
        return 1;
    }

    const auto ref = examples::read_fastx(argv[1]);
    const auto reads = examples::read_fastx(argv[2]);

    // The reference's k-mers with their orientations, from a thread per part
    // of its pieces: `fwd` is whether the k-mer occurs as is.
    const auto pieces = examples::pieces(ref, k);
    kache_hash::Kmer_Map<k, uint8_t> flags(std::max<std::size_t>(examples::positions(ref, k), 1));
    const std::size_t part = (pieces.size() + thread_c - 1) / thread_c;
    std::vector<std::thread> threads;
    for(std::size_t t = 0; t < thread_c; ++t)
        threads.emplace_back([&, t]()
        {
            const auto token = flags.register_user();
            for(std::size_t i = t * part; i < std::min(pieces.size(), (t + 1) * part); ++i)
            {
                const char* const seq = pieces[i].data();
                const std::size_t len = pieces[i].size();
                flags.upsert_oriented(&seq, &len, 1, [](const uint8_t f, const bool fwd){ return uint8_t(f | (fwd ? as_is : revcomp)); }, token);
            }
        });
    for(auto& t : threads)
        t.join();

    // A read's k-mer occurs in the reference's forward strand if it occurs in
    // the reference in its orientation in the read: as is in both, or as the
    // reverse complement in both; and in the reverse strand if in the other
    // one. Each thread looks up a part of the reads, in one batch.
    std::vector<std::size_t> forward(reads.size()), reverse(reads.size());
    const std::size_t read_part = (reads.size() + thread_c - 1) / thread_c;
    threads.clear();
    for(std::size_t t = 0; t < thread_c; ++t)
        threads.emplace_back([&, t]()
        {
            const std::size_t b = std::min(reads.size(), t * read_part), e = std::min(reads.size(), (t + 1) * read_part);
            std::vector<const char*> seqs;
            std::vector<std::size_t> lens;
            for(std::size_t r = b; r < e; ++r)
                seqs.push_back(reads[r].seq.data()), lens.push_back(reads[r].seq.size());

            flags.find(seqs.data(), lens.data(), seqs.size(), [&](const std::size_t q, std::size_t, const std::optional<uint8_t>& f, const bool fwd)
            {
                if(f)
                    forward[b + q] += (*f & (fwd ? as_is : revcomp)) != 0,
                    reverse[b + q] += (*f & (fwd ? revcomp : as_is)) != 0;
            });
        });
    for(auto& t : threads)
        t.join();

    std::cout << "read\tpositions\tforward\treverse\n";
    for(std::size_t r = 0; r < reads.size(); ++r)
    {
        const std::size_t len = reads[r].seq.size();
        std::cout << r << "\t" << (len >= k ? len - k + 1 : 0) << "\t" << forward[r] << "\t" << reverse[r] << "\n";
    }

    return 0;
}
