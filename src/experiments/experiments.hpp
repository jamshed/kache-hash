
#include "Streaming_Kmer_Hash_Table.hpp"
#include "Parser.hpp"

#include "icebergHT/iceberg_table.h"
#include "libcuckoo/cuckoohash_map.hh"
#include "boost/unordered/concurrent_flat_map.hpp"
#include "parlay/parallel.h"
#include "utility.hpp"

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>
#include <cstdlib>
#include <thread>


using namespace kache_hash;


template <uint16_t k, uint16_t l>
uint64_t throughput_kache_hash_sequential_insert(const std::string& path, const std::size_t init_cap)
{
    Streaming_Kmer_Hash_Table<k, false, uint64_t, l> M(init_cap, 1);

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

        Kmer_Window<k, l> w;
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
    std::cerr << "Throughput: " << (n / (duration(t_e - t_s) * 1000'000.0)) << "M op/s.\n";

    return M.size();
}


template <uint16_t k, uint16_t l>
uint64_t throughput_kache_hash_concurrent_insert(const std::string& path, const std::size_t init_cap)
{
    typedef Streaming_Kmer_Hash_Table<k, true, uint64_t, l> ht_t;
    ht_t M(init_cap, parlay::num_workers());

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

            Kmer_Window<k, l> w;
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
    std::cerr << "Throughput: " << (n / (duration(t_e - t_s) * 1000'000.0)) << "M op/s.\n";

    return M.size();
}


template <uint16_t k>
uint64_t throughput_iceberg_concurrent_insert(const std::string& path, const std::size_t init_cap)
{
    iceberg_table M;
    iceberg_init(&M, std::ceil(std::log2(init_cap)));

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
    std::cerr << "Throughput: " << (n / (duration(t_e - t_s) * 1000'000.0)) << "M op/s.\n";

    return tot_balls(&M);
}


template <uint16_t k>
uint64_t throughput_libcuckoo_concurrent_insert(const std::string& path, const std::size_t init_cap)
{
    libcuckoo::cuckoohash_map<Kmer<k>, uint64_t, Kmer_Hasher<k>> M(init_cap);

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
                M.upsert(v.canonical(), [](uint64_t& c)
                {
                    c++;
                }, 1);
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
    std::cerr << "\n Time taken for libcuckoo for " << n << " insertions: " << duration(t_e - t_s) << "s.\n";
    std::cerr << "Throughput: " << (n / (duration(t_e - t_s) * 1000'000.0)) << "M op/s.\n";

    return M.size();
}


template <uint16_t k>
uint64_t throughput_boost_concurrent_insert(const std::string& path, const std::size_t init_cap)
{
    boost::unordered::concurrent_flat_map<Kmer<k>, uint64_t, Kmer_Hasher<k>> M(init_cap);

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
                M.emplace_or_visit(v.canonical(), 1, [](auto& kv)
                {
                    kv.second++;
                });
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
    std::cerr << "\n Time taken for boost for " << n << " insertions: " << duration(t_e - t_s) << "s.\n";
    std::cerr << "Throughput: " << (n / (duration(t_e - t_s) * 1000'000.0)) << "M op/s.\n";

    return M.size();
}


template <uint16_t k, uint16_t l>
void throughput_kache_hash_concurrent_query(const std::string& index_path, const std::string& query_path, const std::size_t init_cap, const std::size_t builder_threads, const uint64_t runs)
{
    typedef Streaming_Kmer_Hash_Table<k, true, uint64_t, l> ht_t;
    ht_t M(init_cap, builder_threads);

{
    uint64_t n = 0;
#ifndef QUIET
    constexpr uint64_t log_step = 32768;
    uint64_t checkpoint = log_step;
    std::cerr << "\n";
#endif

    std::vector<typename ht_t::Token> token(builder_threads);
    for(uint16_t w = 0; w < builder_threads; ++w)
        token[w] = M.register_user();

    Parser input(index_path);
    while(input.read_next_seq())
    {
        std::vector<std::thread> W;
        W.reserve(builder_threads);

        auto const seq = input.seq();
        auto const len = input.seq_len();
        if(len < k)
            continue;

        const auto kmer_per_w = (len - k + 1) / builder_threads;
        for(std::size_t p_id = 0; p_id < builder_threads; ++p_id)
            W.emplace_back([&](const auto p_id)
            {
                const auto s = seq + p_id * kmer_per_w;
                const auto s_len = (p_id < builder_threads - 1 ? kmer_per_w + k - 1 : len - p_id * kmer_per_w);
                if(s_len < k)
                    return;

                Kmer_Window<k, l> w;
                w.init(s);
                std::size_t kmer_idx = 0;
                while(true)
                {
                    M.upsert(w, [](const uint64_t v)
                    {
                        return v + 1;
                    }, 1, token[p_id]);
                    if(kmer_idx + k == s_len)
                        break;

                    w.advance(s[kmer_idx + k]);
                    kmer_idx++;
                }
            }, p_id);

        std::for_each(W.begin(), W.end(), [](auto& w){ w.join(); });

        n += len - k + 1;
#ifndef QUIET
        if(n >= checkpoint)
        {
            std::cerr << "\rInserted " << n << " k-mers.";
            checkpoint = n + log_step;
        }
#endif
    }
}


    for(uint64_t r = 0; r < runs; ++r)
    {
        std::atomic_uint64_t c = 0;
        uint64_t n = 0;
        double t = 0;
    #ifndef QUIET
        constexpr uint64_t log_step = 32768;
        uint64_t checkpoint = log_step;
        std::cerr << "\n";
    #endif

        Parser input(query_path);
        while(input.read_next_seq())
        {
            auto const seq = input.seq();
            auto const len = input.seq_len();
            if(len < k)
                continue;

            const auto t_s = now();
            const auto kmer_per_w = (len - k + 1) / parlay::num_workers();
            parlay::parallel_for(0, parlay::num_workers(), [&](const auto p_id)
            {
                const auto s = seq + p_id * kmer_per_w;
                const auto s_len = (p_id < parlay::num_workers() - 1 ? kmer_per_w + k - 1 : len - p_id * kmer_per_w);
                if(s_len < k)
                    return;

                Kmer_Window<k, l> w;
                w.init(s);
                std::size_t kmer_idx = 0, found = 0;
                while(true)
                {
                    const auto r = M.find(w); assert(r);
                    found += (r ? 1 : 0);

                    if(kmer_idx + k == s_len)
                        break;

                    w.advance(s[kmer_idx + k]);
                    kmer_idx++;
                }

                c += found;
            }
            , 1);

            const auto t_e = now();
            t += duration(t_e - t_s);

            n += len - k + 1;
    #ifndef QUIET
            if(n >= checkpoint)
            {
                std::cerr << "\rQueried " << n << " k-mers.";
                checkpoint = n + log_step;
            }
    #endif
        }

        input.close();
        std::cerr << "\n Time taken for kache-hash for " << c << " positive queries with " << parlay::num_workers() << " threads: " << t << "s at " << " run " << r << ".\n";
        std::cerr << "Throughput: " << (n / (t * 1000'000.0)) << "M op/s.\n";

    }
}


template <uint16_t k>
void throughput_iceberg_concurrent_query(const std::string& index_path, const std::string& query_path, const std::size_t init_cap, const std::size_t builder_threads, const uint64_t runs)
{
    iceberg_table M;
    iceberg_init(&M, std::ceil(std::log2(init_cap)));

{
    uint64_t n = 0;
#ifndef QUIET
    constexpr uint64_t log_step = 32768;
    uint64_t checkpoint = log_step;
    std::cerr << "\n";
#endif

    Parser input(index_path);
    while(input.read_next_seq())
    {
        std::vector<std::thread> W;
        W.reserve(builder_threads);

        auto const seq = input.seq();
        auto const len = input.seq_len();
        if(len < k)
            continue;

        const auto kmer_per_w = (len - k + 1) / builder_threads;
        for(std::size_t p_id = 0; p_id < builder_threads; ++p_id)
            W.emplace_back([&](const auto p_id)
            {
                const auto s = seq + p_id * kmer_per_w;
                const auto s_len = (p_id < builder_threads - 1 ? kmer_per_w + k - 1 : len - p_id * kmer_per_w);
                if(s_len < k)
                    return;

                Directed_Vertex<k> v(s);
                std::size_t kmer_idx = 0;
                while(true)
                {
                    iceberg_insert(&M, v.canonical().as_int(), 1, p_id);
                    if(kmer_idx + k == s_len)
                        break;

                    v.roll_forward(DNA_Utility::map_base(s[kmer_idx + k]));
                    kmer_idx++;
                }
            }
            , p_id);

        std::for_each(W.begin(), W.end(), [](auto& w){ w.join(); });

        n += len - k + 1;

#ifndef QUIET
        if(n >= checkpoint)
        {
            std::cerr << "\rInserted " << n << " k-mers.";
            checkpoint = n + log_step;
        }
#endif
    }
}

    for(uint64_t r = 0; r < runs; ++r)
    {
        std::atomic_uint64_t c = 0;
        uint64_t n = 0;
        double t = 0;
    #ifndef QUIET
        constexpr uint64_t log_step = 32768;
        uint64_t checkpoint = log_step;
        std::cerr << "\n";
    #endif

        Parser input(query_path);
        while(input.read_next_seq())
        {
            std::vector<std::thread> W;
            W.reserve(builder_threads);

            auto const seq = input.seq();
            auto const len = input.seq_len();
            if(len < k)
                continue;

            const auto t_s = now();
            const auto kmer_per_w = (len - k + 1) / parlay::num_workers();
            parlay::parallel_for(0, parlay::num_workers(), [&](const auto p_id)
            {
                const auto s = seq + p_id * kmer_per_w;
                const auto s_len = (p_id < parlay::num_workers() - 1 ? kmer_per_w + k - 1 : len - p_id * kmer_per_w);
                if(s_len < k)
                    return;

                Directed_Vertex<k> v(s);
                std::size_t kmer_idx = 0;
                uint64_t val;
                uint64_t found = 0;
                while(true)
                {
                    val = 0;
                    iceberg_get_value(&M, v.canonical().as_int(), &val, parlay::worker_id());
                    if(val)
                        found++;

                    if(kmer_idx + k == s_len)
                        break;

                    v.roll_forward(DNA_Utility::map_base(s[kmer_idx + k]));
                    kmer_idx++;
                }

                c += found;
            }
            , 1);

            const auto t_e = now();
            t += duration(t_e - t_s);

            n += len - k + 1;
    #ifndef QUIET
            if(n >= checkpoint)
            {
                std::cerr << "\rQueried " << n << " k-mers.";
                checkpoint = n + log_step;
            }
    #endif
        }

        input.close();
        std::cerr << "\n Time taken for iceberg for " << c << " positive queries with " << parlay::num_workers() << " threads: " << t << "s at run " << r << ".\n";
        std::cerr << "Throughput: " << (n / (t * 1000'000.0)) << "M op/s.\n";
    }
}


template <uint16_t k>
void throughput_libcuckoo_concurrent_query(const std::string& index_path, const std::string& query_path, const std::size_t init_cap, const std::size_t builder_threads, const uint64_t runs)
{
    libcuckoo::cuckoohash_map<Kmer<k>, uint64_t, Kmer_Hasher<k>> M(init_cap);

{
    uint64_t n = 0;
#ifndef QUIET
    constexpr uint64_t log_step = 32768;
    uint64_t checkpoint = log_step;
    std::cerr << "\n";
#endif

    Parser input(index_path);
    while(input.read_next_seq())
    {
        std::vector<std::thread> W;
        W.reserve(builder_threads);

        auto const seq = input.seq();
        auto const len = input.seq_len();
        if(len < k)
            continue;

        const auto kmer_per_w = (len - k + 1) / builder_threads;
        for(std::size_t p_id = 0; p_id < builder_threads; ++p_id)
            W.emplace_back([&](const auto p_id)
            {
                const auto s = seq + p_id * kmer_per_w;
                const auto s_len = (p_id < builder_threads - 1 ? kmer_per_w + k - 1 : len - p_id * kmer_per_w);
                if(s_len < k)
                    return;

                Directed_Vertex<k> v(s);
                std::size_t kmer_idx = 0;
                while(true)
                {
                    M.upsert(v.canonical(), [](uint64_t& c)
                    {
                        c++;
                    }, 1);
                    if(kmer_idx + k == s_len)
                        break;

                    v.roll_forward(DNA_Utility::map_base(s[kmer_idx + k]));
                    kmer_idx++;
                }
            }
            , p_id);

        std::for_each(W.begin(), W.end(), [](auto& w){ w.join(); });

        n += len - k + 1;

#ifndef QUIET
        if(n >= checkpoint)
        {
            std::cerr << "\rInserted " << n << " k-mers.";
            checkpoint = n + log_step;
        }
#endif
    }
}


    for(uint64_t r = 0; r < runs; ++r)
    {
        std::atomic_uint64_t c = 0;
        uint64_t n = 0;
        double t = 0;
#ifndef QUIET
        constexpr uint64_t log_step = 32768;
        uint64_t checkpoint = log_step;
        std::cerr << "\n";
#endif

        Parser input(query_path);
        while(input.read_next_seq())
        {
            std::vector<std::thread> W;
            W.reserve(builder_threads);

            auto const seq = input.seq();
            auto const len = input.seq_len();
            if(len < k)
                continue;

            const auto t_s = now();
            const auto kmer_per_w = (len - k + 1) / parlay::num_workers();
            parlay::parallel_for(0, parlay::num_workers(), [&](const auto p_id)
            {
                const auto s = seq + p_id * kmer_per_w;
                const auto s_len = (p_id < parlay::num_workers() - 1 ? kmer_per_w + k - 1 : len - p_id * kmer_per_w);
                if(s_len < k)
                    return;

                Directed_Vertex<k> v(s);
                std::size_t kmer_idx = 0;
                uint64_t found = 0;
                uint64_t val;
                while(true)
                {
                    if(M.find(v.canonical(), val))
                        found++;

                    if(kmer_idx + k == s_len)
                        break;

                    v.roll_forward(DNA_Utility::map_base(s[kmer_idx + k]));
                    kmer_idx++;
                }

                c += found;
            }
            , 1);

            const auto t_e = now();
            t += duration(t_e - t_s);

            n += len - k + 1;
    #ifndef QUIET
            if(n >= checkpoint)
            {
                std::cerr << "\rQueried " << n << " k-mers.";
                checkpoint = n + log_step;
            }
    #endif
        }

        input.close();
        std::cerr << "\n Time taken for libcuckoo for " << c << " positive queries with " << parlay::num_workers() << " threads: " << t << "s at run " << r << ".\n";
        std::cerr << "Throughput: " << (n / (t * 1000'000.0)) << "M op/s.\n";
    }
}


template <uint16_t k>
void throughput_boost_concurrent_query(const std::string& index_path, const std::string& query_path, const std::size_t init_cap, const std::size_t builder_threads, const uint64_t runs)
{
    boost::unordered::concurrent_flat_map<Kmer<k>, uint64_t, Kmer_Hasher<k>> M(init_cap);

{
    uint64_t n = 0;
#ifndef QUIET
    constexpr uint64_t log_step = 32768;
    uint64_t checkpoint = log_step;
    std::cerr << "\n";
#endif

    Parser input(index_path);
    while(input.read_next_seq())
    {
        std::vector<std::thread> W;
        W.reserve(builder_threads);

        auto const seq = input.seq();
        auto const len = input.seq_len();
        if(len < k)
            continue;

        const auto kmer_per_w = (len - k + 1) / builder_threads;
        for(std::size_t p_id = 0; p_id < builder_threads; ++p_id)
            W.emplace_back([&](const auto p_id)
            {
                const auto s = seq + p_id * kmer_per_w;
                const auto s_len = (p_id < builder_threads - 1 ? kmer_per_w + k - 1 : len - p_id * kmer_per_w);
                if(s_len < k)
                    return;

                Directed_Vertex<k> v(s);
                std::size_t kmer_idx = 0;
                while(true)
                {
                    M.emplace_or_visit(v.canonical(), 1, [](auto& kv)
                    {
                        kv.second++;
                    });
                    if(kmer_idx + k == s_len)
                        break;

                    v.roll_forward(DNA_Utility::map_base(s[kmer_idx + k]));
                    kmer_idx++;
                }
            }
            , p_id);

        std::for_each(W.begin(), W.end(), [](auto& w){ w.join(); });

        n += len - k + 1;

#ifndef QUIET
        if(n >= checkpoint)
        {
            std::cerr << "\rInserted " << n << " k-mers.";
            checkpoint = n + log_step;
        }
#endif
    }
}


    for(uint64_t r = 0; r < runs; ++r)
    {
        std::atomic_uint64_t c = 0;
        uint64_t n = 0;
        double t = 0;
    #ifndef QUIET
        constexpr uint64_t log_step = 32768;
        uint64_t checkpoint = log_step;
        std::cerr << "\n";
    #endif

        Parser input(query_path);
        while(input.read_next_seq())
        {
            std::vector<std::thread> W;
            W.reserve(builder_threads);

            auto const seq = input.seq();
            auto const len = input.seq_len();
            if(len < k)
                continue;

            const auto t_s = now();
            const auto kmer_per_w = (len - k + 1) / parlay::num_workers();
            parlay::parallel_for(0, parlay::num_workers(), [&](const auto p_id)
            {
                const auto s = seq + p_id * kmer_per_w;
                const auto s_len = (p_id < parlay::num_workers() - 1 ? kmer_per_w + k - 1 : len - p_id * kmer_per_w);
                if(s_len < k)
                    return;

                Directed_Vertex<k> v(s);
                std::size_t kmer_idx = 0;
                uint64_t found = 0;
                uint64_t val;
                while(true)
                {
                    if(M.visit(v.canonical(), [&val](const auto& kv){ val = kv.second; }))
                        found++;

                    if(kmer_idx + k == s_len)
                        break;

                    v.roll_forward(DNA_Utility::map_base(s[kmer_idx + k]));
                    kmer_idx++;
                }

                c += found;
            }
            , 1);

            const auto t_e = now();
            t += duration(t_e - t_s);

            n += len - k + 1;
    #ifndef QUIET
            if(n >= checkpoint)
            {
                std::cerr << "\rQueried " << n << " k-mers.";
                checkpoint = n + log_step;
            }
    #endif
        }

        input.close();
        std::cerr << "\n Time taken for boost for " << c << " positive queries with " << parlay::num_workers() << " threads: " << t << "s at run " << r << ".\n";
        std::cerr << "Throughput: " << (n / (t * 1000'000.0)) << "M op/s.\n";
    }
}


template <uint16_t k, uint16_t l>
void cache_miss_kache_hash_sequential(const std::string& file_path, const std::size_t init_cap)
{
    throughput_kache_hash_sequential_insert<k, l>(file_path, init_cap);
}


template <uint16_t k, uint16_t l>
void cache_miss_kache_hash_concurrent(const std::string& file_path, const std::size_t init_cap)
{
    throughput_kache_hash_concurrent_insert<k, l>(file_path, init_cap);
}


template <uint16_t k>
void cache_miss_iceberg_concurrent(const std::string& file_path, const std::size_t init_cap)
{
    throughput_iceberg_concurrent_insert<k>(file_path, init_cap);
}


template <uint16_t k>
void cache_miss_libcuckoo_concurrent(const std::string& file_path, const std::size_t init_cap)
{
    throughput_libcuckoo_concurrent_insert<k>(file_path, init_cap);
}


template <uint16_t k>
void cache_miss_boost_concurrent(const std::string& file_path, const std::size_t init_cap)
{
    throughput_boost_concurrent_insert<k>(file_path, init_cap);
}
