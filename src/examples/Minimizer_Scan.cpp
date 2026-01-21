
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
    constexpr uint16_t l = 12;
    Minimizer_Iterator<const char*, k, true> min_it_base(l);
    Min_Iterator<k, l> min_it;

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
                uint64_t min_base, pos_base, h_base;
                min_it_base.value_at(min_base, pos_base, h_base);
                pos_base -= kmer_idx;

                if(min_it.hash() != h_base)
                {
                    std::cerr << "Incorrect minimizer(s) computed. Aborting.\n";
                    std::exit(EXIT_FAILURE);
                }

                if(min_it.hash_pos() >= 0 && min_it.hash_pos() != static_cast<int64_t>(pos_base))
                {
                    std::cerr << "Incorrect position of minimizer(s) in case of unambiguous minimizers. Aborting.\n";
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
    constexpr uint16_t l = 12;
    Min_Iterator<k, l> min_it;

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


void BM_scan_minimizer_pos(benchmark::State& state)
{
    constexpr uint16_t l = 12;
    Min_Iterator<k, l> min_it;

    for(auto _ : state)
    {
        uint64_t min_h = 0;
        int16_t pos_h = 0;
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
                benchmark::DoNotOptimize(pos_h ^= min_it.hash_pos());
                // min_h ^= min_it.hash();
                // pos_h ^= min_it.hash_pos();

                if(kmer_idx + k == len)
                    break;

                min_it.advance(seq[kmer_idx + k]);
                kmer_idx++;
            }
        }

        input.close();

        // std::cout << "Minimizers' signature: " << min_h << "\n";
        // std::cout << "Positions' signature:  " << pos_h << "\n";
    }
}