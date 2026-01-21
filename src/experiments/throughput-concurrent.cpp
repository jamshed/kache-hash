
#include "experiments.hpp"


int main(int argc, char* argv[])
{
    (void)argc;

    constexpr uint16_t k = 31, l = 19;

    std::string path(argv[1]);
    uint64_t init_cap(std::atoll(argv[2]));

    if(argc == 3)
    {
        throughput_kache_hash_concurrent_insert<k, l>(path, init_cap);
        throughput_iceberg_concurrent_insert<k>(path, init_cap);
        throughput_libcuckoo_concurrent_insert<k>(path, init_cap);
        throughput_boost_concurrent_insert<k>(path, init_cap);
    }
    else
    {
        throughput_kache_hash_concurrent_insert_fastq<k, l>(path, init_cap);
        throughput_iceberg_concurrent_insert_fastq<k>(path, init_cap);
        throughput_libcuckoo_concurrent_insert_fastq<k>(path, init_cap);
        throughput_boost_concurrent_insert_fastq<k>(path, init_cap);
    }

    return 0;
}
