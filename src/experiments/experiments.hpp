
#include "Streaming_Kmer_Hash_Table.hpp"
#include "Parser.hpp"

#include "parlay/parallel.h"
#include "utility.hpp"

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>
#include <cstdlib>


using namespace kache_hash;


template <uint16_t k, uint16_t l>
uint64_t throughput_kache_hash_concurrent_insert(const std::string& path, const std::size_t init_cap)
{
    typedef Streaming_Kmer_Hash_Table<k, true, uint64_t, l> ht_t;
    ht_t M(init_cap, parlay::num_workers());

    const auto t_s = now();

    uint64_t n = 0;
#ifndef QUIET
    constexpr uint64_t log_step = 32768;
    uint64_t checkpoint = log_step;
    std::cerr << "\n";
#endif

    std::vector<typename ht_t::Token> token(parlay::num_workers());
    for(uint16_t w = 0; w < parlay::num_workers(); ++w)
        token[w] = M.register_user();

    Parser input(path);
    while(input.read_next_seq())
    {
        auto const seq = input.seq();
        auto const len = input.seq_len();
        if(len < k)
            continue;

        const auto kmer_per_w = (len - k + 1) / parlay::num_workers();
        parlay::parallel_for(0, parlay::num_workers(), [&](const auto p_id)
        {
            const auto s = seq + p_id * kmer_per_w;
            const auto s_len = (p_id < parlay::num_workers() - 1 ? kmer_per_w + k - 1 : len - p_id * kmer_per_w);
            if(s_len < k)
                return;

            Kmer_Window<k, l> w;
            w.init(s);
            std::size_t kmer_idx = 0;
            while(true)
            {
                M.upsert(w, [](const uint64_t v)
                {
                    return v + 1;
                }, 1, token[parlay::worker_id()]);
                if(kmer_idx + k == s_len)
                    break;

                w.advance(s[kmer_idx + k]);
                kmer_idx++;
            }
        }
        , 1);

        n += len - k + 1;
#ifndef QUIET
        if(n >= checkpoint)
        {
            std::cerr << "\rInserted " << n << " k-mers.";
            checkpoint = n + log_step;
        }
#endif
    }

    input.close();
    const auto t_e = now();
    std::cerr << "\n Time taken for kache-hash-concurrent for " << n << " insertions: " << duration(t_e - t_s) << "s.\n";
    std::cerr << "Throughput: " << (n / (duration(t_e - t_s) * 1000'000.0)) << "M op/s.\n";

    return M.size();
}


template <uint16_t k, uint16_t l>
void key_distribution_kache_hash(const std::string& path, const std::size_t init_cap, const double lf = -1)
{
    typedef Streaming_Kmer_Hash_Table<k, true, uint64_t, l> ht_t;
    ht_t M(init_cap, parlay::num_workers(), lf < 0 ? ht_t::lf_default() : lf);
    std::cerr << "Capacity: " << M.capacity() << "\n";

    const auto t_s = now();

    uint64_t n = 0;
#ifndef QUIET
    constexpr uint64_t log_step = 32768;
    uint64_t checkpoint = log_step;
    std::cerr << "\n";
#endif

    std::vector<typename ht_t::Token> token(parlay::num_workers());
    for(uint16_t w = 0; w < parlay::num_workers(); ++w)
        token[w] = M.register_user();

    Parser input(path);
    while(input.read_next_seq())
    {
        auto const seq = input.seq();
        auto const len = input.seq_len();
        if(len < k)
            continue;

        const auto kmer_per_w = (len - k + 1) / parlay::num_workers();
        parlay::parallel_for(0, parlay::num_workers(), [&](const auto p_id)
        {
            const auto s = seq + p_id * kmer_per_w;
            const auto s_len = (p_id < parlay::num_workers() - 1 ? kmer_per_w + k - 1 : len - p_id * kmer_per_w);
            if(s_len < k)
                return;

            Kmer_Window<k, l> w;
            w.init(s);
            std::size_t kmer_idx = 0;
            while(true)
            {
                M.upsert(w, [](const uint64_t v)
                {
                    return v + 1;
                }, 1, token[parlay::worker_id()]);
                if(kmer_idx + k == s_len)
                    break;

                w.advance(s[kmer_idx + k]);
                kmer_idx++;
            }
        }
        , 1);

        n += len - k + 1;
#ifndef QUIET
        if(n >= checkpoint)
        {
            std::cerr << "\rInserted " << n << " k-mers.";
            checkpoint = n + log_step;
        }
#endif
    }

    input.close();
    const auto t_e = now();
    std::cerr << "\n Time taken for kache-hash-concurrent for " << n << " insertions: " << duration(t_e - t_s) << "s.\n";
    std::cerr << "Throughput: " << (n / (duration(t_e - t_s) * 1000'000.0)) << "M op/s.\n";

    std::cerr << "Map size: " << M.size() << "\n";

    const auto d = M.key_distribution();
    std::cerr << "\n\nPercentage of keys going to:\n";
    std::cerr << "\tPrimary:   " << std::get<0>(d) << "%\n";
    std::cerr << "\tSecondary: " << std::get<1>(d) << "%\n";
    std::cerr << "\tOverflow:  " << std::get<2>(d) << "%\n";
}
