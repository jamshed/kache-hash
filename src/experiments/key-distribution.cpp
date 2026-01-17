#include "experiments.hpp"


int main(int argc, char* argv[])
{
    (void)argc;

    constexpr uint16_t k = 31, l = 19;

    std::string path(argv[1]);
    uint64_t init_cap(std::atoll(argv[2]));

    key_distribution_kache_hash<k, l>(path, init_cap);

    return 0;
}
