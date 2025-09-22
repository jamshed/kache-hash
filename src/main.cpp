
#include "Streaming_Kmer_Hashset.hpp"
#include "Parser.hpp"

#include <cstdint>
#include <string>
#include <iostream>
#include <cstdlib>


template <uint16_t k, uint16_t l>
void benchmark_hash_set(const std::string& file_path, const uint64_t final_sz)
{
    Streaming_Kmer_Hashset<k> S(final_sz);

    constexpr uint64_t log_step = 32768;
    uint64_t n = 0, checkpoint = log_step;
    std::cerr << "\n";

    S.clear();

    // std::vector<std::string> P; P.push_back(file_path);
    // S.preprocess_skew_mins(P);
    // S.preprocess_skew_mins_random(P);

    Parser input(file_path);
    while(input.read_next_seq())
    {
        auto const seq = input.seq();
        auto const len = input.seq_len();

        S.template insert<false>(seq, len);

        n += len - k + 1;
        if(n >= checkpoint)
        {
            std::cerr << "\rInserted " << n << " k-mers.";
            checkpoint = n + log_step;
        }
    }

    input.close();
    std::cout << "\nSet size: " << S.size() << "\n";

    std::cerr   << "Complete set size: " << S.size() << ".\n";
    std::cerr   << "Overflow set size: " << S.overflow_size() << " ("
                    << (S.overflow_size() * 100.0 / (S.size() + S.overflow_size())) << "%).\n";
}

int main(int argc, char* argv[])
{
    if(argc != 3)
    {
        std::cerr << "Incorrect number of arguments. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    constexpr uint16_t k = 31, l = 19;
    std::string path(argv[1]);
    uint64_t final_sz(std::atoll(argv[2]));

    std::cerr << "Constructing hash table for sequences from " << path << " with final size " << final_sz << ".\n";

    benchmark_hash_set<k, l>(path, final_sz);

    return 0;
}
