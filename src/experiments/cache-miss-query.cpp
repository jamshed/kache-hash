
#include "experiments.hpp"


int main(int argc, char* argv[])
{
    (void)argc;

    constexpr uint16_t k = 31, l = 17;

    std::string idx_path(argv[1]);
    std::string query_path(argv[2]);
    uint64_t init_cap(std::atoll(argv[3]));
    uint64_t build_threads(std::atoll(argv[4]));
    int opt(std::atoi(argv[5]));

    switch(opt)
    {
    // case 0:
    //     cache_miss_kache_hash_sequential<k, l>(path, init_cap);
    //     break;

    case 1:
        std::cerr << "\n\n======== Baseline run for kache-hash-concurrent ========\n\n";
        throughput_kache_hash_concurrent_query<k, l>(idx_path, query_path, init_cap, build_threads, 1);
        break;

    case 5:
        std::cerr << "\n\n======== Measurement runs for kache-hash-concurrent ========\n\n";
        throughput_kache_hash_concurrent_query<k, l>(idx_path, query_path, init_cap, build_threads, 5);
        break;

    case 2:
        std::cerr << "\n\n======== Baseline run for iceberg ========\n\n";
        throughput_iceberg_concurrent_query<k>(idx_path, query_path, init_cap, build_threads, 1);
        break;

    case 6:
        std::cerr << "\n\n======== Measurement runs for iceberg ========\n\n";
        throughput_iceberg_concurrent_query<k>(idx_path, query_path, init_cap, build_threads, 5);
        break;

    case 3:
        std::cerr << "\n\n======== Baseline run for libcuckoo ========\n\n";
        throughput_libcuckoo_concurrent_query<k>(idx_path, query_path, init_cap, build_threads, 1);
        break;

    case 7:
        std::cerr << "\n\n======== Measurement runs for libcuckoo ========\n\n";
        throughput_libcuckoo_concurrent_query<k>(idx_path, query_path, init_cap, build_threads, 5);
        break;
// /*
    case 4:
        std::cerr << "\n\n======== Baseline run for boost ========\n\n";
        throughput_boost_concurrent_query<k>(idx_path, query_path, init_cap, build_threads, 1);
        break;

    case 8:
        std::cerr << "\n\n======== Measurement runs for boost ========\n\n";
        throughput_boost_concurrent_query<k>(idx_path, query_path, init_cap, build_threads, 5);
        break;
// */

    default:
        std::cerr << "Wrong option.\n";
    }

    return 0;
}

