
// Counts the occurrences in reads of a reference's 31-mers, as PanGenie counts
// its graph's k-mers with Jellyfish, through kache-hash's counter library:
// primes the set with the reference's k-mers, counts the reads' k-mers in it,
// and prints the set's size, the occurrences counted, and the number of the
// set's k-mers by their counts, from 1 to 10. The counter reads its files
// itself, in parallel, and takes pipes too.
// Usage: count_in_set <reference.fa> <reads.fa|reads.fq> [threads]

#include <kache_hash/Kmer_Set_Counter.hpp>

#include <cstddef>
#include <filesystem>
#include <iostream>
#include <string>


constexpr uint16_t k = 31;
constexpr std::size_t max_count = 10;


int main(const int argc, const char* const argv[])
{
    if(argc < 3)
    {
        std::cerr << "usage: count_in_set <reference.fa> <reads.fa|reads.fq> [threads]\n";
        return 1;
    }

    const std::size_t thread_c = argc > 3 ? std::stoul(argv[3]) : 4;
    if(thread_c == 0 || thread_c > 128)
    {
        std::cerr << "threads: 1 to 128\n";
        return 1;
    }

    // The reference's bytes bound its k-mers; the counter grows beyond them as
    // needed, as for a pipe's.
    std::error_code ec;
    const auto bytes = std::filesystem::file_size(argv[1], ec);
    kache_hash::Kmer_Set_Counter counter(k, ec ? std::size_t(1) << 22 : bytes, thread_c);

    counter.prime(argv[1]);
    const auto counted = counter.count(argv[2]);
    const auto hist = counter.histogram(max_count);

    std::cout << "set " << k << "-mers\t" << counter.size() << "\n";
    std::cout << "occurrences counted\t" << counted << "\n";
    for(std::size_t c = 1; c <= max_count; ++c)
        std::cout << "count " << c << "\t" << hist[c] << "\n";

    return 0;
}
