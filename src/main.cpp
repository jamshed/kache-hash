
#include "DNA_Utility.hpp"
#include "Directed_Vertex.hpp"
#include "Rolling_Hash.hpp"
#include "Streaming_Kmer_Hashset.hpp"
#include "Parser.hpp"
#include "globals.hpp"

#include "nthash/nthash.hpp"
#include "unordered_dense/unordered_dense.h"
#include "parlay/parallel.h"

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iostream>


template <uint16_t k, uint16_t l>
void benchmark_hash_set(const std::string& file_path, const uint64_t final_sz)
{
    Streaming_Kmer_Hash_Table<k, l> S(final_sz);

    constexpr uint64_t log_step = 32768;
    uint64_t n = 0, checkpoint = log_step;
    std::cerr << "\n";

    S.clear();

    Parser input(file_path);
    while(input.read_next_seq())
    {
        auto const seq = input.seq();
        auto const len = input.seq_len();

        S.template insert<false>(seq, len);

        n += len - k + 1;
        if(n >= checkpoint)
        {
            std::cerr << "\rInserted " << n << " k-mers.";
            checkpoint = n + log_step;
        }
    }

    input.close();
    std::cout << "\nSet size: " << S.size() << "\n";

    std::cerr   << "Complete set size: " << S.size() << ".\n";
    std::cerr   << "Overflow set size: " << S.overflow_size() << " ("
                    << (S.overflow_size() * 100.0 / (S.size() + S.overflow_size())) << "%).\n";
}


template <uint16_t k, bool canonical>
void test_rolling_hash(const std::string& file_path)
{
    std::vector<uint64_t> C(256, 0);
    ankerl::unordered_dense::map<Kmer<k>, uint64_t, Kmer_Hasher<k>> kmer_hash;
    ankerl::unordered_dense::set<uint64_t> H;

    uint64_t kmer_c = 0;
    Parser input(file_path);
    while(input.read_next_seq())
    {
        auto const seq = input.seq();
        auto const len = input.seq_len();
        assert(len >= k);

        Directed_Vertex<k> v(seq);
        Rolling_Hash<k, canonical> rh(seq);
        std::size_t kmer_idx = 0;
        while(true)
        {
            const Kmer<k>& kmer = canonical ? v.canonical() : v.kmer();
            if(!kmer_hash.contains(kmer))
                kmer_hash.emplace(kmer, rh.hash()),
                H.insert(rh.hash()),
                C[rh.checksum()]++;
            else if(kmer_hash[kmer] != rh.hash())
            {
                std::cerr << "Hash mismatch for k-mer " << kmer << " at index " << kmer_idx << ". Previously computed hash: " << kmer_hash[kmer] << ", newly computed hash: " << rh.hash() << "\n";
                std::exit(EXIT_FAILURE);
            }

            kmer_c++;
            if(kmer_idx + k >= len)
                break;

            v.roll_forward(DNA_Utility::map_base(seq[kmer_idx + k]));
            rh.advance(seq[kmer_idx + k]);
            kmer_idx++;
        }
    }

    std::cerr << "\n======== Hash stats ========\n";
    std::cerr << "#Inserted-k-mers: " << kmer_c << "\n";
    std::cerr << "#Unique-k-mers: " << kmer_hash.size() << "\n";
    std::cerr << "#Unique-hashes: " << H.size() << "\n";
    std::cerr << "Hash uniqueness ratio: " << ((H.size() * 100.0) / kmer_hash.size()) << "%.\n";

    std::cerr << "\n======== Checksum stats ========\n";
    double c_p2 = 0;
    uint64_t max_x = 0, min_x = std::numeric_limits<uint64_t>::max(), sum_x = 0;
    std::for_each(C.cbegin(), C.cend(),
        [&](auto x) { c_p2 += x * x, max_x = std::max(max_x, x), min_x = std::min(min_x, x), sum_x += x; });

    const auto exp = sum_x / 256.0;
    const double var = (c_p2 / 256.0 - exp * exp);
    std::cerr << "Avg-load: " << static_cast<uint64_t>(exp + 0.5) << "\n";
    std::cerr << "Max-load: " << max_x << "\n";
    std::cerr << "Min-load: " << min_x << "\n";
    std::cerr << "SD(load): " << std::sqrt(var) << "\n";
}


template <uint16_t k>
void check_rolling_hash_reverse(const std::string& file_path)
{
    std::vector<uint64_t> C(256, 0);
    ankerl::unordered_dense::map<Kmer<k>, uint64_t, Kmer_Hasher<k>> kmer_hash;
    ankerl::unordered_dense::set<uint64_t> H;

    Parser input(file_path);
    while(input.read_next_seq())
    {
        auto const seq = input.seq();
        auto const len = input.seq_len();

        std::string rc;
        rc.reserve(len);
        for(std::size_t i = 0; i < len; ++i)
            rc.push_back(DNA_Utility::complement(seq[len - 1 - i]));

        std::vector<uint64_t> H_f;
        std::vector<uint64_t> H_r;

        {
            Rolling_Hash<k, true> rh(seq);
            std::size_t i = 0;
            while(true)
            {
                H_f.push_back(rh.hash_rev());
                if(i + k >= len)
                    break;

                rh.advance(seq[i + k]);
                i++;
            }
        }

        {
            Rolling_Hash<k, false> rh_r(rc.data());
            std::size_t i = 0;
            while(true)
            {
                H_r.push_back(rh_r.hash());
                if(i + k >= len)
                    break;

                rh_r.advance(rc[i + k]);
                i++;
            }
        }

        std::reverse(H_r.begin(), H_r.end());
        for(std::size_t i = 0; i < H_f.size(); ++i)
            if(H_f[i] != H_r[i])
            {
                std::cerr << "Hashes do not match for index " << i << ".\n";
                std::cerr << "Forward hash: " << H_f[i] << ", reverse hash: " << H_r[i] << "\n";
                return;
            }
    }
}


template <uint16_t k>
void bench_rolling_hash(const std::string& file_path)
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
            h ^= rh.hash();

            if(kmer_idx + k == len)
                break;

            rh.advance(seq[kmer_idx + k]);
            kmer_idx++;
        }
    }

    std::cout << "Input k-mer collection signature: " << h << "\n";

    input.close();
}


template <uint16_t k>
void bench_ntHash(const std::string& file_path)
{
    Parser input(file_path);
    uint64_t h = 0;

    while(input.read_next_seq())
    {
        auto const seq = input.seq();
        // auto const len = input.seq_len();
        // assert(len >= k);
        // Rolling_Hash<k, true> rh(seq);
        nthash::NtHash nth(seq, 1, k);

        // std::size_t kmer_idx = 0;
        while(nth.roll())
        {
            h ^= *nth.hashes();

            // if(kmer_idx + k == len)
            //     break;

            // rh.advance(seq[kmer_idx + k]);
            // kmer_idx++;
        }
    }

    std::cout << "Input k-mer collection signature: " << h << "\n";

    input.close();
}


template <uint16_t k, uint16_t l>
void test_mt_insert()
{
    Streaming_Kmer_Hash_Table<k, l> S(final_sz);

    std::size_t n = 0;
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

        const auto kmer_per_w = (len - (k - 1)) / parlay::num_workers();
        parlay::parallel_for(0, parlay::num_workers(), [&](const auto w_id)
        {
            const auto s = seq + w_id * kmer_per_w;
            const auto s_len = (w_id < parlay::num_workers() - 1 ? kmer_per_w + k - 1 : len - w_id * kmer_per_w);
            Kmer_Window<k> w;

            w.init(s);
            std::size_t kmer_idx = 0;
            while(true)
            {
                S.insert_mt(w);
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

    std::cout << "\nSet size: " << S.template size<true>() << "\n";
    std::cout << "Overflow-set size: " << S.template overflow_size<true>() << "\n";
}


int main(int argc, char* argv[])
{
    (void)argc, (void)argv;

    constexpr uint16_t k = 31, l = 19;
    // (void)k, (void)l;

    // std::string path(argv[1]);
    // uint64_t final_sz(std::atoll(argv[2]));

    // benchmark_hash_set<k, l>(path, final_sz);

    // test_rolling_hash<k, true>(path);
    // test_rolling_hash<k, true>(path);
    // check_rolling_hash_reverse<k>(path);
    // bench_rolling_hash<k>(path);
    // bench_ntHash<k>(path);

    test_mt_insert<k, l>();

    return 0;
}
