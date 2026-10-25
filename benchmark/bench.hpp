
#ifndef KACHE_HASH_BENCH_HPP
#define KACHE_HASH_BENCH_HPP



#include "harness.hpp"
#include "utility.hpp"

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <thread>
#include <algorithm>
#include <optional>
#include <memory>


// A table to benchmark is driven through an adapter type `Table_` that has:
// - `name`, `is_set`, and `concurrent`: static constants for its name, whether
//   it is a set rather than a map of k-mer counts, and whether it is safe for
//   multiple user threads;
// - `report_params()`: a static function printing its parameters;
// - `Table_(init_cap, thread_c)`: a constructor sizing it for `init_cap` k-mers
//   and `thread_c` user threads;
// - `insert(t_id, s, n)`: adds the `n` k-mers of the sequence `s`, for the
//   `t_id`'th user thread—incrementing their counts in a map;
// - `find(t_id, s, n)`: returns the number of the `n` k-mers of the sequence `s`
//   present in it, for the `t_id`'th user thread;
// - `size()`: its number of k-mers;
// - `capacity()`: its number of k-mer slots—only if it has slots;
// - `report()`: prints its statistics.


// Benchmark configuration, from the command line.
struct Config
{
    std::string index;  // Path to the input to build the table on.
    std::vector<std::string> query; // Paths to the inputs to query the table with.
    std::size_t thread_c = 1;   // Number of user threads.
    bool seq = false;   // Whether to use the sequential variant of the table.
    bool set = false;   // Whether to use a set instead of a map of k-mer counts.
    std::size_t init_cap = 0;   // Number of k-mers to initially size the table for; the index k-mer count if 0.
    std::optional<std::string> perf;    // Raw hardware events to count per phase, as `name=rCODE` items, if any.
};


// Returns the benchmark configuration from the command line arguments `argv`,
// exiting on invalid ones.
inline Config parse_config(const int argc, char* argv[])
{
    const auto usage = [&]()
    {
        std::cerr << "Usage: " << argv[0] << " --index <path> [--query <path>]... [--threads <count>] [--seq] [--set] [--init-cap <count>] [--perf <events>]\n";
        std::exit(EXIT_FAILURE);
    };

    Config cfg;
    for(int i = 1; i < argc; ++i)
    {
        const std::string arg(argv[i]);
        const bool has_val = (i + 1 < argc);
        if(arg == "--index" && has_val)
            cfg.index = argv[++i];
        else if(arg == "--query" && has_val)
            cfg.query.emplace_back(argv[++i]);
        else if(arg == "--threads" && has_val)
            cfg.thread_c = std::strtoull(argv[++i], nullptr, 10);
        else if(arg == "--seq")
            cfg.seq = true;
        else if(arg == "--set")
            cfg.set = true;
        else if(arg == "--init-cap" && has_val)
            cfg.init_cap = std::strtoull(argv[++i], nullptr, 10);
        else if(arg == "--perf" && has_val)
            cfg.perf = argv[++i];
        else
            usage();
    }

    if(cfg.index.empty() || cfg.thread_c == 0 || cfg.thread_c > std::thread::hardware_concurrency() || (cfg.seq && cfg.thread_c > 1))
        usage();

    return cfg;
}


// Exits on the configuration being unsupported by the table, for the reason
// `reason`.
inline void reject(const std::string& reason)
{
    std::cerr << "Unsupported configuration: " << reason << ".\n";
    std::exit(EXIT_FAILURE);
}


// Applies `f` to the canonical form of each of the `n` k-mers of the sequence
// `s`.
template <typename F_>
void for_each_canonical(const char* const s, const std::size_t n, const F_& f)
{
    Directed_Vertex<k> v(Kmer<k>{s});
    std::size_t kmer_idx = 0;
    while(true)
    {
        f(v.canonical());
        if(kmer_idx + 1 == n)
            break;

        v.roll_forward(DNA_Utility::map_base(s[kmer_idx + k]));
        kmer_idx++;
    }
}


// Benchmarks the table `Table_` per the configuration `cfg`: builds it over the
// k-mers of the index, and then looks up the k-mers of each query input in it.
// Counts the configured hardware events per phase, if any.
template <typename Table_>
void bench(const Config& cfg)
{
    // The counters only follow the threads spawned after they are opened.
    const auto perf = (cfg.perf ? std::make_unique<Perf_Group>(*cfg.perf) : nullptr);

    print("structure", Table_::name);
    print("k", k);
    Table_::report_params();
    print("table", Table_::is_set ? "set" : "map");
    print("concurrent", Table_::concurrent);
    print("threads", cfg.thread_c);
    if(cfg.perf)
        print("perf.events", *cfg.perf);

    const auto t_load_s = now();
    const Fragments index(cfg.index);
    const auto t_load_e = now();
    report_input("index", cfg.index, index);
    print("load.seconds", duration(t_load_e - t_load_s));
    report_rss("input");

    const auto init_cap = (cfg.init_cap > 0 ? cfg.init_cap : index.kmer_count());
    const auto t_init_s = now();
    Table_ M(init_cap, cfg.thread_c);
    const auto t_init_e = now();
    print("table.init_cap", init_cap);
    if constexpr(requires { M.capacity(); })
        print("table.init_capacity", M.capacity());
    print("init.seconds", duration(t_init_e - t_init_s));

    const auto t_build = run_parallel(index, cfg.thread_c, [&](const std::size_t t_id, const char* const s, const std::size_t n)
    {
        M.insert(t_id, s, n);
    }, perf.get());

    report_phase("build", index.kmer_count(), t_build);
    if(perf)
        report_perf("build", *perf);

    report_rss("build");
    if constexpr(requires { M.capacity(); })
        print("table.capacity", M.capacity());
    print("table.size", M.size());
    M.report();

    for(std::size_t q = 0; q < cfg.query.size(); ++q)
    {
        const auto name = "query" + std::to_string(q);
        const Fragments query(cfg.query[q]);
        report_input(name, cfg.query[q], query);

        std::vector<Padded<std::size_t>> H(cfg.thread_c, Padded<std::size_t>(0));
        const auto t_query = run_parallel(query, cfg.thread_c, [&](const std::size_t t_id, const char* const s, const std::size_t n)
        {
            H[t_id].unwrap() += M.find(t_id, s, n);
        }, perf.get());

        std::size_t hit_c = 0;
        std::for_each(H.cbegin(), H.cend(), [&](const auto& h){ hit_c += h.unwrap(); });
        report_phase(name, query.kmer_count(), t_query);
        print(name + ".hits", hit_c);
        if(perf)
            report_perf(name, *perf);
    }

    report_memory();
}



#endif
