
#include "Parser.hpp"
#include "Kmer.hpp"
#include "Directed_Vertex.hpp"
#include "Rolling_Hash.hpp"
#include "Streaming_Kmer_Hash_Table.hpp"
#include "utility.hpp"

// #include "nthash/nthash.hpp"
#include "unordered_dense/unordered_dense.h"

#include <cstdint>
#include <cstdlib>
#include <iostream>


uint64_t benchmark_file_read(const std::string& path)
{
    uint64_t sum_len = 0;
    Parser input(path);
    while(input.read_next_seq())
        sum_len += input.seq_len();

    input.close();
    return sum_len;
}


uint64_t benchmark_file_scan(const std::string& path)
{
    uint64_t a = 0, c = 0, g = 0, t = 0;
    Parser input(path);
    while(input.read_next_seq())
    {
        const auto s = input.seq();
        for(std::size_t i = 0; i < input.seq_len(); ++i)
            a += (s[i] == 'A'),
            c += (s[i] == 'C'),
            g += (s[i] == 'G'),
            t += (s[i] == 'T');
    }

    input.close();
    return a ^ c ^ g ^ t;
}


template <uint16_t k>
uint64_t benchmark_kmer_hash(const std::string& path)
{
    Parser input(path);
    uint64_t h = 0;

    while(input.read_next_seq())
    {
        auto const seq = input.seq();
        auto const len = input.seq_len();
        if(len < k)
            continue;

        Directed_Vertex<k> v(seq);
        std::size_t kmer_idx = 0;
        while(true)
        {
            h ^= v.canonical().to_u64();

            if(kmer_idx + k == len)
                break;

            v.roll_forward(DNA_Utility::map_base(seq[kmer_idx + k]));
            kmer_idx++;
        }
    }

    input.close();
    return h;
}


template <uint16_t k>
uint64_t benchmark_rolling_hash(const std::string& path)
{
    Parser input(path);
    uint64_t h = 0;

    while(input.read_next_seq())
    {
        auto const seq = input.seq();
        auto const len = input.seq_len();
        if(len < k)
            continue;

        Rolling_Hash<k, true> rh(seq);
        std::size_t kmer_idx = 0;
        while(true)
        {
            h ^= rh.hash();

            if(kmer_idx + k == len)
                break;

            rh.advance(seq[kmer_idx + k]);
            kmer_idx++;
        }
    }

    input.close();
    return h;
}


/*
template <uint16_t k>
uint64_t benchmark_ntHash(const std::string& path)
{
    Parser input(path);
    uint64_t h = 0;

    while(input.read_next_seq())
    {
        auto const seq = input.seq();
        if(input.seq_len() < k)
            continue;

        nthash::NtHash nth(seq, 1, k);

        while(nth.roll())
            h ^= *nth.hashes();
    }

    input.close();
    return h;
}
*/


template <uint16_t k, typename T_ht_>
uint64_t benchmark_set_insert_generic_sequential_no_resize(const std::string& path, const std::size_t final_sz)
{
    T_ht_ S(final_sz);

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
    std::cerr << "\n Time: " << duration(t_e - t_s) << "s.\n";

    return S.size();
}


template <uint16_t k>
uint64_t benchmark_set_insert_streaming_sequential_no_resize(const std::string& path, const std::size_t final_sz)
{
    Streaming_Kmer_Hash_Table<k, false> S(final_sz);

    const auto t_s = now();

    uint64_t n = 0;
#ifndef QUIET
    constexpr uint64_t log_step = 32768;
    uint64_t checkpoint = log_step;
    std::cerr << "\n";
#endif

    Parser input(path);
    Kmer_Window<k> w;
    while(input.read_next_seq())
    {
        auto const seq = input.seq();
        auto const len = input.seq_len();
        if(len < k)
            continue;

        w.init(seq);
        std::size_t kmer_idx = 0;
        while(true)
        {
            S.insert(w);
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
    std::cerr << "\n Time: " << duration(t_e - t_s) << "s.\n";

    return S.size();
}


int main(int argc, char* argv[])
{
    (void)argc, (void)argv;
    if(argc != 3)
    {
        std::cerr << "Pass two arguments: file-path and final size of hash table.\n";
        return EXIT_FAILURE;
    }

    constexpr uint16_t k = 31, l = 19;

    std::string path(argv[1]);
    uint64_t final_sz(std::atoll(argv[2]));
    (void)l, (void)final_sz;

    // std::cerr << benchmark_file_read(path) << "\n";
    // std::cerr << benchmark_file_scan(path) << "\n";
    // std::cerr << benchmark_kmer_hash<k>(path) << "\n";
    // std::cerr << benchmark_rolling_hash<k>(path) << "\n";
    // std::cerr << benchmark_ntHash<k>(path) << "\n";
    // std::cerr << benchmark_set_insert_generic_sequential_no_resize<k, ankerl::unordered_dense::set<Kmer<k>, Kmer_Hasher<k>>>(path, final_sz) << "\n";
    std::cerr << benchmark_set_insert_streaming_sequential_no_resize<k>(path, final_sz) << "\n";


    return EXIT_SUCCESS;
}
