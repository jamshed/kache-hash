
#include "experiments.hpp"


int main(int argc, char* argv[])
{
    (void)argc;

    constexpr uint16_t k = 31, l = 17;
    constexpr uint64_t runs = 4;

    std::string idx_path(argv[1]);
    std::string query_path(argv[2]);
    uint64_t init_cap(std::atoll(argv[3]));
    uint64_t build_threads(std::atoll(argv[4]));

    throughput_kache_hash_concurrent_query<k, l>(idx_path, query_path, init_cap, build_threads, runs);
    throughput_iceberg_concurrent_query<k>(idx_path, query_path, init_cap, build_threads, runs);
    throughput_libcuckoo_concurrent_query<k>(idx_path, query_path, init_cap, build_threads, runs);
    throughput_boost_concurrent_query<k>(idx_path, query_path, init_cap, build_threads, runs);

    return 0;
}

