// A consumer of the counting library: the 31-mers of a sequence, counted within the set of its own.

#include <kache_hash/Kmer_Set_Counter.hpp>

#include <cstdio>
#include <fstream>
#include <string>


int main()
{
    const std::string path = "consumer-counter.fa";
    std::ofstream(path) << ">s\nATATCACACCCAACCTTCAAATGCCGTGCCCTAACGCCCTAATCCTGCGC\n";  // 20 distinct 31-mers.

    kache_hash::Kmer_Set_Counter counter(31, 1000, 2);
    counter.prime(path);
    const auto counted = counter.count(path);

    std::printf("kache_counter: %zu 31-mers in the set, %zu occurrences counted\n", counter.size(), counted);
    return counter.size() == 20 && counted == 20 ? 0 : 1;
}
