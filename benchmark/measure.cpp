
#include "harness.hpp"
#include "Streaming_Kmer_Hash_Table.hpp"
#include "Op_Stats.hpp"
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
    bool scan = false;  // Whether to scan the index k-mers without the table.
    std::size_t thread_c = 1;   // Number of user threads.
    bool seq = false;   // Whether to use the non-concurrent table.
    bool set = false;   // Whether to use a hash set instead of a hash map of k-mer counts.
    std::size_t init_cap = 0;   // Number of k-mers to initially size the table for; the index k-mer count if 0.
    std::optional<std::string> perf;    // Raw hardware events to count per phase, as `name=rCODE` items, if any.
};


// Scans the k-mers in the fragments `F` with `thread_c` threads, computing all
// that the table operations extract from the k-mer windows, without the table:
// the canonical k-mer, its checksum, and its minimizer hash and coordinate.
// Puts a digest of these at `digest`, and returns the elapsed time in seconds.
// Counts the events of the group `perf` over the scan, if provided.
double scan(const Fragments& F, const std::size_t thread_c, std::size_t& digest, Perf_Group* const perf)
{
    std::vector<Padded<std::size_t>> D(thread_c, Padded<std::size_t>(0));
    const auto t = run_parallel(F, thread_c, [&](const std::size_t t_id, const char* const s, const std::size_t n)
    {
        // Mirrors `Kmer_Window`, so that none of its computations are optimized out.
        Directed_Vertex<k> v(Kmer<k>{s});
        Rolling_Hash<k, true> rh;
        Min_Iterator<k, l> min_it;
        rh.init(s), min_it.reset(s);

        std::size_t d = 0;
        std::size_t kmer_idx = 0;
        while(true)
        {
            uint8_t m;
            const auto c = std::max(rh.checksum<8>(), uint64_t(1));
            const auto h = min_it.hash(m);
            m = v.in_canonical_form() ? m : (m ^ Min_Iterator<k, l>::min_orientation_mask());
            d += v.canonical().as_int() ^ h ^ (c << 8) ^ m;
            if(kmer_idx + 1 == n)
                break;

            const auto ch = s[kmer_idx + k];
            v.roll_forward(DNA_Utility::map_base(ch));
            min_it.advance(ch);
            rh.advance(ch);
            kmer_idx++;
        }

        D[t_id].unwrap() += d;
    }, perf);

    digest = 0;
    std::for_each(D.cbegin(), D.cend(), [&](const auto& d){ digest += d.unwrap(); });
    return t;
}


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
    if constexpr(instrument)
        print("instrument", instrument);
    if constexpr(no_table_lock)
        print("no_table_lock", no_table_lock);

    const Fragments index(cfg.index);
    report_input("index", cfg.index, index);

    if(cfg.scan)
    {
        std::size_t digest;
        const auto t = scan(index, cfg.thread_c, digest, perf);
        report_phase("scan", index.kmer_count(), t);
        print("scan.digest", digest);
        if(perf)
            report_perf("scan", *perf);
    }

    const auto init_cap = (cfg.init_cap > 0 ? cfg.init_cap : index.kmer_count());
    ht_t M(init_cap, cfg.thread_c);
    print("table.init_cap", init_cap);
    print("table.init_capacity", M.capacity());

    std::vector<typename ht_t::Token> token(cfg.thread_c);
    std::for_each(token.begin(), token.end(), [&](auto& t){ t = M.register_user(); });

    if constexpr(instrument)
        Op_Stats::reset();

    const auto t_build = run_parallel(index, cfg.thread_c, [&](const std::size_t t_id, const char* const s, const std::size_t n)
    {
        for_each_window(s, n, [&](const Kmer_Window<k, l>& w)
        {
            if constexpr(is_set)
                M.insert(w, token[t_id]);
            else
                M.upsert(w, [](const uint64_t v){ return v + 1; }, 1, token[t_id]);
        });
    }, perf);

    report_phase("build", index.kmer_count(), t_build);
    if(perf)
        report_perf("build", *perf);
    if constexpr(instrument)
        Op_Stats::total().report(std::cout, "inst.build.");

    print("table.capacity", M.capacity());
    print("table.size", M.size());
    print("table.overflow_size", M.overflow_size());

    const auto [primary, secondary, overflow] = M.key_distribution();
    print("table.primary_pct", primary);
    print("table.secondary_pct", secondary);
    print("table.overflow_pct", overflow);
    if constexpr(instrument)
    {
        print("table.overflow_capacity", M.overflow_capacity());
        const auto occ = M.occupancy();
        for(std::size_t s = 0; s < occ.size(); ++s)
            print("table.occupancy." + std::to_string(s), occ[s]);
    }

    for(std::size_t q = 0; q < cfg.query.size(); ++q)
    {
        const auto name = "query" + std::to_string(q);
        const Fragments query(cfg.query[q]);
        report_input(name, cfg.query[q], query);

        std::vector<Padded<std::size_t>> H(cfg.thread_c, Padded<std::size_t>(0));
        if constexpr(instrument)
            Op_Stats::reset();

        const auto t_query = run_parallel(query, cfg.thread_c, [&](const std::size_t t_id, const char* const s, const std::size_t n)
        {
            std::size_t h = 0;
            for_each_window(s, n, [&](const Kmer_Window<k, l>& w){ h += static_cast<bool>(M.find(w)); });
            H[t_id].unwrap() += h;
        }, perf);

        std::size_t hit_c = 0;
        std::for_each(H.cbegin(), H.cend(), [&](const auto& h){ hit_c += h.unwrap(); });
        report_phase(name, query.kmer_count(), t_query);
        print(name + ".hits", hit_c);
        if(perf)
            report_perf(name, *perf);
        if constexpr(instrument)
            Op_Stats::total().report(std::cout, "inst." + name + ".");
    }

    report_memory();
}


// Returns the measurement configuration from the command line arguments
// `argv`, exiting on invalid ones.
Config parse_config(const int argc, char* argv[])
{
    const auto usage = [&]()
    {
        std::cerr << "Usage: " << argv[0] << " --index <path> [--query <path>]... [--scan] [--threads <count>] [--seq] [--set] [--init-cap <count>] [--perf <events>]\n";
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
        else if(arg == "--scan")
            cfg.scan = true;
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

    // The table supports up to 64 users, and no more than the hardware threads.
    const std::size_t max_thread_c = std::min(64u, std::thread::hardware_concurrency());
    if(cfg.index.empty() || cfg.thread_c == 0 || cfg.thread_c > max_thread_c || (cfg.seq && cfg.thread_c > 1))
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
