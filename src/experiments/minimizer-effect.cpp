#include "experiments.hpp"


int main(int argc, char* argv[])
{
    (void)argc;

    constexpr uint16_t k = 31;

    const std::string path(argv[1]);
    const uint64_t init_cap(std::atoll(argv[2]));

#ifdef MIN_SIZE_EXP
    throughput_kache_hash_concurrent_insert<k, MIN_SIZE_EXP>(path, init_cap);
#else
    (void)k, (void)init_cap;
    std::cerr << "`MIN_SIZE_EXP` is undefined.\n";
#endif

    return 0;
}
