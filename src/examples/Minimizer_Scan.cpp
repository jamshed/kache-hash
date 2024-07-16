
#include "Minimizer_Iterator.hpp"
#include "Parser.hpp"
#include "globals.hpp"
#include "benchmark/benchmark.h"

#include <cstdint>
#include <cstddef>
#include <cassert>


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
