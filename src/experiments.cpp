
#include "Streaming_Kmer_Hash_Table.hpp"
#include "Parser.hpp"

#include "icebergHT/iceberg_table.h"
#include "parlay/parallel.h"
#include "utility.hpp"

#include <cstdint>
#include <cstddef>
#include <string>

template <uint16_t k>
void throughput_kache_hash_sequential(const std::string& path, const std::size_t final_sz)
{
    Streaming_Kmer_Hash_Table<k, false, uint64_t> M(final_sz, 1);

    const auto t_s = now();

    uint64_t n = 0;
#ifndef QUIET
    constexpr uint64_t log_step = 32768;
    uint64_t checkpoint = log_step;
    std::cerr << "\n";
#endif

    Parser input(path);
    const auto token = M.register_user();
    while(input.read_next_seq())
    {
        auto const seq = input.seq();
        auto const len = input.seq_len();
        if(len < k)
            continue;

        Kmer_Window<k> w;
        w.init(seq);
        std::size_t kmer_idx = 0;
        while(true)
        {
            M.upsert(w, [](const uint64_t v)
            {
                return v + 1;
            }, 1, token);

            if(kmer_idx + k == len)
                break;

            w.advance(seq[kmer_idx + k]);
            kmer_idx++;
        }

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
    std::cerr << "\n Time taken for kache-hash-sequential for " << n << " insertions: " << duration(t_e - t_s) << "s.\n";
    std::cerr << "Throughput: " << duration(t_e - t_s) / (n * 1000'000.0) << "M op/s.\n";

    return M.size();
}


template <uint16_t k>
uint64_t throughput_kache_hash_concurrent(const std::string& path, const std::size_t final_sz)
{
    typedef Streaming_Kmer_Hash_Table<k, true, uint64_t> ht_t;
    ht_t M(final_sz, parlay::num_workers());

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

            Kmer_Window<k> w;
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
    std::cerr << "Throughput: " << duration(t_e - t_s) / (n * 1000'000.0) << "M op/s.\n";

    return M.size();
}


template <uint16_t k>
uint64_t benchmark_iceberg_map_insert_generic_concurrent_no_resize(const std::string& path, const std::size_t final_sz)
{
    iceberg_table M;
    iceberg_init(&M, std::ceil(std::log2(final_sz)));

    const auto t_s = now();

    uint64_t n = 0;
#ifndef QUIET
    constexpr uint64_t log_step = 32768;
    uint64_t checkpoint = log_step;
    std::cerr << "\n";
#endif

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

            Directed_Vertex<k> v(s);
            std::size_t kmer_idx = 0;
            while(true)
            {
                iceberg_insert(&M, v.canonical().as_int(), 1, parlay::worker_id());
                if(kmer_idx + k == s_len)
                    break;

                v.roll_forward(DNA_Utility::map_base(s[kmer_idx + k]));
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
    std::cerr << "\n Time taken for iceberg for " << n << " insertions: " << duration(t_e - t_s) << "s.\n";
    std::cerr << "Throughput: " << duration(t_e - t_s) / (n * 1000'000.0) << "M op/s.\n";

    return tot_balls(&M);
}


int main(int argc, char* argv[])
{
    (void)argc, (void)argv;

    return 0;
}
