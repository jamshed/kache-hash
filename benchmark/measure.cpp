
#include "harness.hpp"
#include "Streaming_Kmer_Hash_Table.hpp"
#include "utility.hpp"

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>
#include <type_traits>
#include <iostream>
#include <cstdlib>
#include <thread>
#include <algorithm>
#include <optional>
#include <memory>


using namespace kache_hash;


// Measurement configuration, from the command line.
struct Config
{
    std::string index;  // Path to the input to build the table on.
    std::vector<std::string> query; // Paths to the inputs to query the table with.
    std::size_t thread_c = 1;   // Number of user threads.
    bool seq = false;   // Whether to use the non-concurrent table.
    bool set = false;   // Whether to use a hash set instead of a hash map of k-mer counts.
    std::size_t init_cap = 0;   // Number of k-mers to initially size the table for; the index k-mer count if 0.
    std::optional<std::string> perf;    // Raw hardware events to count per phase, as `name=rCODE` items, if any.
    std::size_t reps = 1;   // Number of repetitions of the build and the queries, each on a fresh table.
};


constexpr std::size_t batch_sz = 64;    // Most pieces of fragments per call of the streaming operations.


// Measures the table operations per the configuration `cfg`, on a table that is
// concurrent iff `mt_`, and that maps k-mers to `T_` values—a set if `T_` is
// `void`. Counts the events of the group `perf` per phase, if provided.
template <bool mt_, typename T_>
void measure(const Config& cfg, Perf_Group* const perf)
{
    typedef Streaming_Kmer_Hash_Table<k, mt_, T_, l> ht_t;
    constexpr bool is_set = std::is_void_v<T_>;

    print("k", k);
    print("l", l);
    print("table", is_set ? "set" : "map");
    print("concurrent", mt_);
    print("threads", cfg.thread_c);
    if(cfg.perf)
        print("perf.events", *cfg.perf);

    const Fragments index(cfg.index);
    report_input("index", cfg.index, index);

    std::vector<std::unique_ptr<const Fragments>> query_own;    // Query inputs other than the index.
    std::vector<const Fragments*> query;    // `query[q]` is the `q`'th query input.
    for(std::size_t q = 0; q < cfg.query.size(); ++q)
    {
        if(cfg.query[q] != cfg.index)
            query_own.emplace_back(std::make_unique<const Fragments>(cfg.query[q]));
        query.push_back(cfg.query[q] == cfg.index ? &index : query_own.back().get());
        report_input("query" + std::to_string(q), cfg.query[q], *query.back());
    }

    for(std::size_t rep = 0; rep < cfg.reps; ++rep)
    {
        print("rep", rep);

        const auto init_cap = (cfg.init_cap > 0 ? cfg.init_cap : index.kmer_count());
        const auto t_init_s = now();
        ht_t M(init_cap, cfg.thread_c);
        const auto t_init_e = now();
        print("table.init_cap", init_cap);
        print("table.init_capacity", M.capacity());
        print("init.seconds", duration(t_init_e - t_init_s));

        std::vector<typename ht_t::Token> token(cfg.thread_c);
        std::for_each(token.begin(), token.end(), [&](auto& t){ t = M.register_user(); });

        const auto t_build = run_parallel_batched(index, cfg.thread_c, batch_sz, [&](const std::size_t t_id, const char* const* const seqs, const std::size_t* const lens, const std::size_t c)
        {
            if constexpr(is_set)
                M.insert(seqs, lens, c, token[t_id]);
            else
                M.upsert(seqs, lens, c, [](const uint64_t v){ return v + 1; }, 1, token[t_id]);
        }, perf);

        report_phase("build", index.kmer_count(), t_build);
        if(perf)
            report_perf("build", *perf);

        print("table.capacity", M.capacity());
        print("table.size", M.size());
        print("table.overflow_size", M.overflow_size());

        const auto [primary, secondary, overflow] = M.key_distribution();
        print("table.primary_pct", primary);
        print("table.secondary_pct", secondary);
        print("table.overflow_pct", overflow);

        for(std::size_t q = 0; q < query.size(); ++q)
        {
            const auto name = "query" + std::to_string(q);
            std::vector<Padded<std::size_t>> H(cfg.thread_c, Padded<std::size_t>(0));
            const auto t_query = run_parallel_batched(*query[q], cfg.thread_c, batch_sz, [&](const std::size_t t_id, const char* const* const seqs, const std::size_t* const lens, const std::size_t c)
            {
                H[t_id].unwrap() += M.count(seqs, lens, c);
            }, perf);

            std::size_t hit_c = 0;
            std::for_each(H.cbegin(), H.cend(), [&](const auto& h){ hit_c += h.unwrap(); });
            report_phase(name, query[q]->kmer_count(), t_query);
            print(name + ".hits", hit_c);
            if(perf)
                report_perf(name, *perf);
        }

        report_memory();    // With the table still in place.
    }
}


// Returns the measurement configuration from the command line arguments
// `argv`, exiting on invalid ones.
Config parse_config(const int argc, char* argv[])
{
    const auto usage = [&]()
    {
        std::cerr << "Usage: " << argv[0] << " --index <path> [--query <path>]... [--threads <count>] [--seq] [--set] [--init-cap <count>] [--perf <events>] [--reps <count>]\n";
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
        else if(arg == "--reps" && has_val)
            cfg.reps = std::strtoull(argv[++i], nullptr, 10);
        else
            usage();
    }

    // The table supports up to `max_user` users, and no more than the hardware threads.
    const std::size_t max_thread_c = std::min<std::size_t>(Streaming_Kmer_Hash_Table<k, true, uint64_t, l>::max_user, std::thread::hardware_concurrency());
    if(cfg.index.empty() || cfg.thread_c == 0 || cfg.reps == 0 || cfg.thread_c > max_thread_c || (cfg.seq && cfg.thread_c > 1))
        usage();

    return cfg;
}


int main(int argc, char* argv[])
{
    const auto cfg = parse_config(argc, argv);

    // The counters only follow the threads spawned after they are opened.
    const auto perf = (cfg.perf ? std::make_unique<Perf_Group>(*cfg.perf) : nullptr);

    if(cfg.set)
        cfg.seq ? measure<false, void>(cfg, perf.get()) : measure<true, void>(cfg, perf.get());
    else
        cfg.seq ? measure<false, uint64_t>(cfg, perf.get()) : measure<true, uint64_t>(cfg, perf.get());

    return 0;
}
