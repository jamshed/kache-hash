
#include "DNA_Utility.hpp"
#include "Kmer.hpp"
#include "Directed_Vertex.hpp"
#include "Rolling_Hash.hpp"
#include "Parser.hpp"
#include "Kmer_Hashset.hpp"
#include "Streaming_Double_Probe_Hashset.hpp"
#include "Streaming_Hashset.hpp"
#include "Streaming_Kmer_Hashset.hpp"
#include "globals.hpp"

#include "nthash/nthash.hpp"
#include "unordered_dense/unordered_dense.h"
#include "benchmark/benchmark.h"
#include "parlay/parallel.h"

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>


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


void BM_file_scan(benchmark::State& state)
{
    uint64_t a, c, g, t;
    for(auto _ : state)
    {
        a = 0, c = 0, g = 0, t = 0;
        Parser input(file_path);
        while(input.read_next_seq())
        {
            const auto s = input.seq();
            for(std::size_t i = 0; i < input.seq_len(); ++i)
                benchmark::DoNotOptimize(a += (s[i] == 'A')),
                benchmark::DoNotOptimize(c += (s[i] == 'C')),
                benchmark::DoNotOptimize(g += (s[i] == 'G')),
                benchmark::DoNotOptimize(t += (s[i] == 'T'));
        }
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


void BM_rolling_hash(benchmark::State& state)
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
            Rolling_Hash<k, true> rh(seq);

            std::size_t kmer_idx = 0;
            while(true)
            {
                benchmark::DoNotOptimize(h ^= rh.hash());

                if(kmer_idx + k == len)
                    break;

                rh.advance(seq[kmer_idx + k]);
                kmer_idx++;
            }
        }

        // std::cout << "Input k-mer set signature: " << h << "\n";

        input.close();
    }
}


void BM_ntHash(benchmark::State& state)
{
    for(auto _ : state)
    {
        Parser input(file_path);
        uint64_t h = 0;

        while(input.read_next_seq())
        {
            auto const seq = input.seq();
            assert(input.seq_len() >= k);
            nthash::NtHash nth(seq, 1, k);

            while(nth.roll())
                benchmark::DoNotOptimize(h ^= *nth.hashes());
        }

        // std::cout << "Input k-mer set signature: " << h << "\n";

        input.close();
    }
}


template <typename T_ht_>
void BM_set_insert_generic(benchmark::State& state)
{
    T_ht_ S(final_sz);
    // ankerl::unordered_dense::set<Kmer<k>, Kmer_Hasher<k>> S;
    // S.reserve(final_sz);

    for(auto _ : state)
    {
        uint64_t n = 0;
#ifndef QUIET
        constexpr uint64_t log_step = 32768;
        uint64_t checkpoint = log_step;
        std::cerr << "\n";
#endif

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
                // S.insert(v.canonical());
                benchmark::DoNotOptimize(S.insert(v.canonical()));

                if(kmer_idx + k == len)
                    break;

                v.roll_forward(DNA_Utility::map_base(seq[kmer_idx + k]));
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

        std::cout << "\nSet size: " << S.size() << "\n";
    }
}


void BM_set_insert_generic_3rd_party(benchmark::State& state)
{
    BM_set_insert_generic<ankerl::unordered_dense::set<Kmer<k>, Kmer_Hasher<k>>>(state);
}


void BM_set_insert_generic(benchmark::State& state)
{
    BM_set_insert_generic<Kmer_Hashset<k>>(state);
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
                // S.insert_lazy(v.canonical());
                benchmark::DoNotOptimize(S.insert_lazy(v.canonical()));

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
    Streaming_Double_Probe_Hashset<k> S(final_sz);

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


template <bool track_stats_>
void BM_set_insert_streaming_bucketed(benchmark::State& state)
{
    Streaming_Kmer_Hashset<k> S(final_sz);

    for(auto _ : state)
    {
        uint64_t n = 0;
#ifndef QUIET
        constexpr uint64_t log_step = 32768;
        uint64_t checkpoint = log_step;
        std::cerr << "\n";
#endif

        S.clear();

        Parser input(file_path);
        while(input.read_next_seq())
        {
            auto const seq = input.seq();
            auto const len = input.seq_len();

            S.template insert<track_stats_>(seq, len);

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
    }

    std::cout << "\nSet size: " << S.size<false>() << "\n";
    std::cout << "Overflow-set size: " << S.overflow_size<false>() << "\n";
    std::ofstream os(ht_path, std::ios::out | std::ios::binary);
    S.serialize(os);
    os.close();

    if constexpr(track_stats_)
        S.report_stats();
}


void BM_set_insert_streaming_bucketed_mt(benchmark::State& state)
{
    Streaming_Kmer_Hashset<k> S(final_sz);

    for(auto _ : state)
    {
        std::atomic_uint64_t n = 0;
#ifndef QUIET
        constexpr uint64_t log_step = 32768;
        uint64_t checkpoint = log_step;
        std::cerr << "\n";
#endif

        S.clear();

        Parser input(file_path);
        while(input.read_next_seq())
        {
            auto const seq = input.seq();
            auto const len = input.seq_len();

            const auto kmer_per_w = len / parlay::num_workers();
            parlay::parallel_for(0, parlay::num_workers(), [&](const auto p_id)
            {
                const auto sz = (p_id < parlay::num_workers() - 1 ? kmer_per_w + k - 1 : len - p_id * kmer_per_w);
                S.insert_mt(seq + p_id * kmer_per_w, sz);
                n += sz - k + 1;
            }
            , 1);

#ifndef QUIET
            if(n >= checkpoint)
            {
                std::cerr << "\rInserted " << n << " k-mers.";
                checkpoint = n + log_step;
            }
#endif
        }

        input.close();
    }

    std::cout << "\nSet size: " << S.size<true>() << "\n";
    std::ofstream os(ht_path, std::ios::out | std::ios::binary);
    S.serialize(os);
    os.close();
}


void BM_set_query_positive_streaming(benchmark::State& state)
{
    std::ifstream is(ht_path, std::ios::in | std::ios::binary);
    Streaming_Kmer_Hashset<k> S(is);

    std::cerr << "Set size: " << S.size<false>() << "\n";

    for(auto _ : state)
    {
        uint64_t c = 0, n = 0;
#ifndef QUIET
        constexpr uint64_t log_step = 32768;
        uint64_t checkpoint = log_step;
        std::cerr << "\n";
#endif

        Parser input(file_path);
        while(input.read_next_seq())
        {
            auto const seq = input.seq();
            auto const len = input.seq_len();

            c += S.find(seq, len);

            n += len - k + 1;
#ifndef QUIET
            if(n >= checkpoint)
            {
                std::cerr << "\rQueried " << n << " k-mers.";
                checkpoint = n + log_step;
            }
#endif
        }

        std::cerr << "Positive query: " << c << ".\n";
    }

    is.close();
}



template void BM_set_insert_streaming_blocked<false>(benchmark::State& state);
template void BM_set_insert_streaming_blocked<true>(benchmark::State& state);
template void BM_set_insert_streaming_bucketed<true>(benchmark::State& state);
template void BM_set_insert_streaming_bucketed<false>(benchmark::State& state);
