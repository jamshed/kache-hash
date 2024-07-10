
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


template <uint16_t k>
inline void set_insert_generic(const std::string& file_path)
{
    ankerl::unordered_dense::set<Kmer<k>, Kmer_Hasher<k>> S;

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


void BM_set_insert_generic(benchmark::State& state)
{
    constexpr uint16_t k = 31;
    const std::string file_path("/ssd2/jamshed/small-read.k31.cf3.fa");
    for(auto _ : state)
        set_insert_generic<k>(file_path);
};
