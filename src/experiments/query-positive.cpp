
#include "experiments.hpp"


int main(int argc, char* argv[])
{
    (void)argc;

    constexpr uint16_t k = 31, l = 17;

    std::string path(argv[1]);
    uint64_t init_cap(std::atoll(argv[2]));
    uint64_t build_threads(std::atoll(argv[3]));

    throughput_kache_hash_concurrent_query<k, l>(path, init_cap, build_threads);
    throughput_iceberg_concurrent_query<k>(path, init_cap, build_threads);
    throughput_libcuckoo_concurrent_query<k>(path, init_cap, build_threads);
    throughput_boost_concurrent_query<k>(path, init_cap, build_threads);

    return 0;
}

