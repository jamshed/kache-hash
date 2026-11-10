// A consumer of the header-only table: a set of 31-mers, through the library's umbrella header.

#include <kache_hash/kache_hash.hpp>

#include <cstdio>
#include <cstring>


int main()
{
    const char* const seq = "ATATCACACCCAACCTTCAAATGCCGTGCCCTAACGCCCTAATCCTGCGC";  // 20 distinct 31-mers.
    const std::size_t len = std::strlen(seq);

    kache_hash::Kmer_Set<31> set(1000);
    const auto token = set.register_user();
    const auto inserted = set.insert(seq, len, token);
    const auto present = set.count(seq, len);

    std::printf("kache-hash %s: %zu 31-mers inserted, %zu present\n", KACHE_HASH_VERSION_STRING, inserted, present);
    return inserted == 20 && present == 20 ? 0 : 1;
}
