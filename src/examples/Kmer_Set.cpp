
#include "DNA_Utility.hpp"
#include "Kmer.hpp"
#include "Directed_Vertex.hpp"
#include "Parser.hpp"
#include "unordered_dense/unordered_dense.h"
#include "benchmark/benchmark.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <iostream>
#include <cassert>


void BM_set_insert_generic(benchmark::State& state)
{
    constexpr uint16_t k = 31;
    const std::string file_path("/ssd2/jamshed/small-read.k31.cf3.fa");
    constexpr std::size_t final_sz = 8501727;

    ankerl::unordered_dense::set<Kmer<k>, Kmer_Hasher<k>> S;
    S.reserve(final_sz);

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
                S.insert(v.canonical());

                if(kmer_idx + k == len)
                    break;

                v.roll_forward(DNA_Utility::map_base(seq[kmer_idx + k]));
                kmer_idx++;
            }
        }

        input.close();

        std::cout << "Set size: " << S.size() << "\n";
    }
};
