
#include "experiments.hpp"


int main(int argc, char* argv[])
{
    (void)argc;

    constexpr uint16_t k = 31, l = 17;

    std::string path(argv[1]);
    uint64_t init_cap(std::atoll(argv[2]));
    int opt(std::atoi(argv[3]));

    switch(opt)
    {
    case 0:
        cache_miss_kache_hash_sequential<k, l>(path, init_cap);
        break;

    case 1:
        cache_miss_kache_hash_concurrent<k, l>(path, init_cap);
        break;

    case 2:
        cache_miss_iceberg_concurrent<k>(path, init_cap);
        break;

    case 3:
        cache_miss_libcuckoo_concurrent<k>(path, init_cap);
        break;

    case 4:
        cache_miss_boost_concurrent<k>(path, init_cap);
        break;

    default:
        std::cerr << "Wrong option.\n";
    }

    return 0;
}

