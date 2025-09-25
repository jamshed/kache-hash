
#include "DNA_Utility.hpp"
#include "Rolling_Hash.hpp"
#include "Streaming_Kmer_Hashset.hpp"
#include "Parser.hpp"

#include "unordered_dense/unordered_dense.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <iostream>
#include <cmath>
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


template <uint16_t k, uint16_t l>
void test_rolling_hash(const std::string& file_path)
{
    std::vector<uint64_t> C(256, 0);
    ankerl::unordered_dense::map<Kmer<k>, uint64_t, Kmer_Hasher<k>> kmer_hash;
    ankerl::unordered_dense::set<uint64_t> H;

    Parser input(file_path);
    while(input.read_next_seq())
    {
        auto const seq = input.seq();
        auto const len = input.seq_len();

        Kmer<k> kmer(seq);
        Rolling_Hash<k> rh(seq);
        std::size_t i = 0;
        while(true)
        {

            if(!kmer_hash.contains(kmer))
                kmer_hash.emplace(kmer, rh.hash()),
                H.insert(rh.hash()),
                C[rh.checksum()]++;
            else if(kmer_hash[kmer] != rh.hash())
            {
                std::cerr << "Hash mismatch for k-mer " << kmer << " at index " << i << ". Previously computed hash: " << kmer_hash[kmer] << ", newly computed hash: " << rh.hash() << "\n";
                std::exit(EXIT_FAILURE);
            }

            if(i + k >= len)
                break;

            kmer.roll_forward(DNA_Utility::map_base(seq[i + k]));
            rh.advance(seq[i + k]);
            i++;
        }
    }

    std::cerr << "\n======== Hash stats ========\n";
    std::cerr << "#Unique-k-mers: " << kmer_hash.size() << "\n";
    std::cerr << "#Unique-hashes: " << H.size() << "\n";
    std::cerr << "Hash uniqueness ratio: " << ((H.size() * 100.0) / kmer_hash.size()) << "%.\n";

    std::cerr << "\n======== Checksum stats ========\n";
    double c_p2 = 0;
    uint64_t max_x = 0, min_x = std::numeric_limits<uint64_t>::max(), sum_x = 0;
    std::for_each(C.cbegin(), C.cend(),
        [&](auto x) { c_p2 += x * x, max_x = std::max(max_x, x), min_x = std::min(min_x, x), sum_x += x; });

    const auto exp = sum_x / 256.0;
    const double var = (c_p2 / 256.0 - exp * exp);
    std::cerr << "Avg-load: " << static_cast<uint64_t>(exp + 0.5) << "\n";
    std::cerr << "Max-load: " << max_x << "\n";
    std::cerr << "Min-load: " << min_x << "\n";
    std::cerr << "SD(load): " << std::sqrt(var) << "\n";
}

int main(int argc, char* argv[])
{
    (void)argc, (void)argv;

    constexpr uint16_t k = 31, l = 19;
    std::string path(argv[1]);
    // uint64_t final_sz(std::atoll(argv[2]));

    // benchmark_hash_set<k, l>(path, final_sz);

    test_rolling_hash<k, l>(path);

    return 0;
}
