
#include "Minimizer_Iterator_Baseline.hpp"
#include "Minimizer_Iterator.hpp"
#include "Parser.hpp"
#include "globals.hpp"
#include "benchmark/benchmark.h"

#include <cstdint>
#include <cstddef>
#include <iostream>
#include <cstdlib>
#include <cassert>


void scan_minimizer_correctness(benchmark::State& state)
{
    const uint16_t l = 12;
    Minimizer_Iterator<const char*, k, true> min_it_base(l);
    Min_Iterator<k> min_it(l);

    for(auto _ :state)
    {
        Parser input(file_path);
        while(input.read_next_seq())
        {
            auto const seq = input.seq();
            auto const len = input.seq_len();
            assert(len >= k);

            std::size_t kmer_idx = 0;
            min_it_base.reset(seq, len);
            min_it.reset(seq);

            while(true)
            {
                if(min_it.hash() != min_it_base.hash())
                {
                    std::cerr << "Incorrect minimizer(s) computed. Aborting.\n";
                    std::exit(EXIT_FAILURE);
                }

                if(kmer_idx + k == len)
                    break;

                min_it_base.advance(seq[kmer_idx + k]);
                min_it.advance(seq[kmer_idx + k]);
                kmer_idx++;
            }
        }

        input.close();
    }


    std::cerr << "Minimizers are correct.\n";
}


void BM_scan_minimizer_baseline(benchmark::State& state)
{
    const uint16_t l = 12;
    Minimizer_Iterator<const char*, k, true> min_it(l);

    for(auto _ : state)
    {
        uint64_t min_h = 0;
        Parser input(file_path);
        while(input.read_next_seq())
        {
            auto const seq = input.seq();
            auto const len = input.seq_len();
            assert(len >= k);

            std::size_t kmer_idx = 0;
            min_it.reset(seq, len);

            while(true)
            {
                benchmark::DoNotOptimize(min_h ^= min_it.hash());
                // min_h ^= min_it.hash();

                if(kmer_idx + k == len)
                    break;

                min_it.advance(seq[kmer_idx + k]);
                kmer_idx++;
            }
        }

        input.close();

        // std::cout << "Minimizers' signature: " << min_h << "\n";
    }
}


void BM_scan_minimizer(benchmark::State& state)
{
    const uint16_t l = 12;
    Min_Iterator<k> min_it(l);

    for(auto _ : state)
    {
        uint64_t min_h = 0;
        Parser input(file_path);
        while(input.read_next_seq())
        {
            auto const seq = input.seq();
            auto const len = input.seq_len();
            assert(len >= k);

            std::size_t kmer_idx = 0;
            min_it.reset(seq);

            while(true)
            {
                benchmark::DoNotOptimize(min_h ^= min_it.hash());
                // min_h ^= min_it.hash();

                if(kmer_idx + k == len)
                    break;

                min_it.advance(seq[kmer_idx + k]);
                kmer_idx++;
            }
        }

        input.close();

        // std::cout << "Minimizers' signature: " << min_h << "\n";
    }
}
