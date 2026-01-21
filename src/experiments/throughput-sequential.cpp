
#include "experiments.hpp"


int main(int argc, char* argv[])
{
    (void)argc;

    constexpr uint16_t k = 31, l = 17;

    std::string path(argv[1]);
    uint64_t init_cap(std::atoll(argv[2]));

    throughput_kache_hash_sequential_insert<k, l>(path, init_cap);
    throughput_iceberg_concurrent_insert<k>(path, init_cap);
    throughput_libcuckoo_concurrent_insert<k>(path, init_cap);
    throughput_boost_concurrent_insert<k>(path, init_cap);

    return 0;
}
