#include "experiments.hpp"


template <uint16_t K, uint16_t L>
void run(const uint16_t l, const std::string &path, const std::size_t init_cap)
{
    if constexpr(L == 0)
        return;
    else if(L == l)
        key_distribution_kache_hash<K, 17>(path, init_cap);
    else if(L >= K)
        run<K, K - 1>(l, path, init_cap);
    else
        run<K, L - 1>(l, path, init_cap);
}


int main(int argc, char* argv[])
{
    (void)argc;

    constexpr uint16_t k = 31;

    const std::string path(argv[1]);
    const uint64_t init_cap(std::atoll(argv[2]));

#ifdef MIN_SIZE_EXP
    throughput_kache_hash_concurrent_insert<k, MIN_SIZE_EXP>(path, init_cap);
#else
    std::cerr << "`MIN_SIZE_EXP` is undefined.\n";
#endif

    return 0;
}
