
#include "DNA_Utility.hpp"
#include "Kmer.hpp"
#include "Directed_Vertex.hpp"
#include "Parser.hpp"
#include "Kmer_Hashset.hpp"
#include "Streaming_Kmer_Hashset.hpp"
#include "Streaming_Hashset.hpp"
#include "globals.hpp"
#include "unordered_dense/unordered_dense.h"
#include "benchmark/benchmark.h"

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <iostream>
#include <cassert>
#include <sys/types.h>


void BM_file_read(benchmark::State& state)
{
    for(auto _ : state)
    {
        std::size_t sum_len = 0;
        Parser input(file_path);
        while(input.read_next_seq())
            benchmark::DoNotOptimize(sum_len += input.seq_len());

        // std::cout << "Sequence length: " << sum_len << "\n";
    }
}


void BM_kmer_hash(benchmark::State& state)
{
    for(auto _ : state)
    {
        Parser input(file_path);
        uint64_t h = 0;

        while(input.read_next_seq())
        {
            auto const seq = input.seq();
            auto const len = input.seq_len();
            assert(len >= k);

            Directed_Vertex<k> v(seq);
            std::size_t kmer_idx = 0;
            while(true)
            {
                benchmark::DoNotOptimize(h ^= v.canonical().to_u64());

                if(kmer_idx + k == len)
                    break;

                v.roll_forward(DNA_Utility::map_base(seq[kmer_idx + k]));
                kmer_idx++;
            }
        }

        // std::cout << "Input k-mer set signature: " << h << "\n";

        input.close();
    }
}


void BM_set_insert_generic(benchmark::State& state)
{
    // ankerl::unordered_dense::set<Kmer<k>, Kmer_Hasher<k>> S;
    // S.reserve(final_sz);
    Kmer_Hashset<k> S(final_sz);

    for(auto _ : state)
    {
        S.clear();

        Parser input(file_path);
        while(input.read_next_seq())
        {
            auto const seq = input.seq();
            auto const len = input.seq_len();
            assert(len >= k);

            Directed_Vertex<k> v(seq);
            std::size_t kmer_idx = 0;
            while(true)
            {
                // TODO: this statement is significantly slower than the next one w/ the custom hashset. Diagnose.
                S.insert(v.canonical());
                // benchmark::DoNotOptimize(S.insert(v.canonical()));

                if(kmer_idx + k == len)
                    break;

                v.roll_forward(DNA_Utility::map_base(seq[kmer_idx + k]));
                kmer_idx++;
            }
        }

        input.close();

        std::cout << "Set size: " << S.size() << "\n";
    }
}


void BM_set_insert_batched(benchmark::State& state)
{
    Kmer_Hashset<k> S(final_sz);

    for(auto _ : state)
    {
        S.clear();

        Parser input(file_path);
        while(input.read_next_seq())
        {
            auto const seq = input.seq();
            auto const len = input.seq_len();
            assert(len >= k);

            Directed_Vertex<k> v(seq);
            std::size_t kmer_idx = 0;
            while(true)
            {
                // TODO: this statement is significantly slower than the next one w/ the custom hashset. Diagnose.
                S.insert_lazy(v.canonical());
                // benchmark::DoNotOptimize(S.insert_lazy(v.canonical()));

                if(kmer_idx + k == len)
                    break;

                v.roll_forward(DNA_Utility::map_base(seq[kmer_idx + k]));
                kmer_idx++;
            }
        }

        S.flush_buffer();
        input.close();

        std::cout << "Set size: " << S.size() << "\n";
    }
}


void BM_set_insert_streaming(benchmark::State& state)
{
    Streaming_Kmer_Hashset<k> S(final_sz);

    for(auto _ : state)
    {
        S.clear();

        Parser input(file_path);
        while(input.read_next_seq())
        {
            auto const seq = input.seq();
            auto const len = input.seq_len();
            assert(len >= k);

            S.insert(seq, len);
        }

        input.close();

        std::cout << "Set size: " << S.size() << "\n";
    }
}


template <bool track_stats_>
void BM_set_insert_streaming_blocked(benchmark::State& state)
{
    Streaming_Hash_Set<k> S(final_sz);

    for(auto _ : state)
    {
        S.clear();

        Parser input(file_path);
        while(input.read_next_seq())
        {
            auto const seq = input.seq();
            auto const len = input.seq_len();
            assert(len >= k);

            S.template insert<track_stats_>(seq, len);
        }

        input.close();
    }

    std::cout << "Set size: " << S.size() << "\n";
    if constexpr(track_stats_)
        S.print_stats();
}



template void BM_set_insert_streaming_blocked<false>(benchmark::State& state);
template void BM_set_insert_streaming_blocked<true>(benchmark::State& state);
