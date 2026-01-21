#include "experiments.hpp"
#include <cstdlib>


int main(int argc, char* argv[])
{
    (void)argc;

    constexpr uint16_t k = 31;
#ifdef MIN_SIZE_EXP
    constexpr uint16_t l = MIN_SIZE_EXP;
#else
    constexpr uint16_t l = 19;
#endif

    const std::string path(argv[1]);
    const uint64_t init_cap(std::atoll(argv[2]));
    const double lf(std::atof(argv[3]));

    key_distribution_kache_hash<k, l>(path, init_cap, lf);

    return 0;
}
