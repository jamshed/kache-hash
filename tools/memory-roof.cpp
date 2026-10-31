#include "Perf_Group.hpp"

#include <cstdint>
#include <cstddef>
#include <cstdlib>
#include <string>
#include <vector>
#include <thread>
#include <atomic>
#include <chrono>
#include <optional>
#include <memory>
#include <algorithm>
#include <utility>
#include <iostream>
#include <fstream>

#include <sys/mman.h>


// Measures the memory roofs of a roofline analysis on a machine, with a given
// number of threads:
//  (a) the rates of sequential reads,
//  (b) of independent random reads and updates of cache lines,
//  (c) and of dependent random reads, using a pointer chase whose rate per
//      thread is the inverse of the memory latency, over an array as large as
//      the structures that the roofs bound.
// Each is a phase of a fixed duration, counting a group of hardware events if
// given. It prints a line per measurement, a key and its value.


constexpr std::size_t line_words = L1_CACHE_LINE_SIZE / sizeof(uint64_t);  // Number of 64-bit words in a cache line.


// Configuration of the roofs' measurement.
struct Config
{
    std::size_t thread_c = 1;   // Number of threads.
    std::size_t bytes = std::size_t(32) << 30;  // Size of the array: a power of two, of at least a page.
    std::size_t depth = 32; // Number of random accesses that each thread prefetches ahead of its next one.
    double seconds = 2; // Duration of each phase.
    std::optional<std::string> perf;    // Raw hardware events to count per phase, as `name=rCODE` items, if any.
};


// A count of one thread, on a cache line of its own.
struct alignas(L1_CACHE_LINE_SIZE) Count
{
    std::size_t val = 0;
};


// Returns the configuration of the roofs' measurement from the command line
// arguments `argv`, exiting on invalid ones.
Config parse_config(const int argc, char* argv[])
{
    const auto usage = [&]()
    {
        std::cerr << "Usage: " << argv[0] << " [--threads <count>] [--bytes <count>] [--depth <count>] [--seconds <seconds>] [--perf <events>]\n";
        std::exit(EXIT_FAILURE);
    };

    Config cfg;
    for(int i = 1; i < argc; ++i)
    {
        const std::string arg(argv[i]);
        const bool has_val = (i + 1 < argc);
        if(arg == "--threads" && has_val)
            cfg.thread_c = std::strtoull(argv[++i], nullptr, 10);
        else if(arg == "--bytes" && has_val)
            cfg.bytes = std::strtoull(argv[++i], nullptr, 10);
        else if(arg == "--depth" && has_val)
            cfg.depth = std::strtoull(argv[++i], nullptr, 10);
        else if(arg == "--seconds" && has_val)
            cfg.seconds = std::strtod(argv[++i], nullptr);
        else if(arg == "--perf" && has_val)
            cfg.perf = argv[++i];
        else
            usage();
    }

    const bool pow2 = (cfg.bytes & (cfg.bytes - 1)) == 0;
    if(cfg.thread_c == 0 || cfg.thread_c > std::thread::hardware_concurrency() || !pow2 || cfg.bytes < 4096 || cfg.depth == 0 || !(cfg.seconds > 0))
        usage();

    return cfg;
}


// Prints the measurement `key` with value `val`.
template <typename T_>
void print(const std::string& key, const T_& val)
{
    std::cout << key << "\t" << val << std::endl;
}


// Prints the measurements of the phase `name` that accessed `line_c` cache
// lines in `t` seconds.
void report_phase(const std::string& name, const std::size_t line_c, const double t)
{
    print(name + ".ops", line_c);
    print(name + ".seconds", t);
    print(name + ".mops", line_c / (t * 1'000'000.0));
    print(name + ".gb_per_s", line_c * L1_CACHE_LINE_SIZE / (t * 1e9));
}


// Prints the counts of the performance-counter events in the group `perf` for
// the phase `name`.
void report_perf(const std::string& name, const Perf_Group& perf)
{
    for(std::size_t i = 0; i < perf.size(); ++i)
        print("perf." + name + "." + perf.name(i), perf.count(i));

    print("perf." + name + ".running_pct", perf.running_pct());
}


// Prints the anonymous memory of the process, and how much of it transparent
// huge pages back, as of the point `name` in the run; and the mode of the
// transparent huge pages.
void report_huge_pages(const std::string& name)
{
#ifdef __linux__
    std::ifstream smaps("/proc/self/smaps_rollup");
    std::string field;
    std::size_t kb;
    while(smaps >> field)
        if(field == "Anonymous:" && smaps >> kb)
            print("memory." + name + "_anon_kb", kb);
        else if(field == "AnonHugePages:" && smaps >> kb)
            print("memory." + name + "_anon_huge_pages_kb", kb);

    std::ifstream thp("/sys/kernel/mm/transparent_hugepage/enabled");   // Lists the modes, e.g. "always [madvise] never".
    std::string mode;
    while(thp >> mode)
        if(mode.front() == '[')
            print("memory.thp", mode.substr(1, mode.size() - 2));
#else
    (void)name;
#endif
}


// Returns the next pseudo-random number of the xorshift* generator of state
// `x`, advancing it.
inline uint64_t next_random(uint64_t& x)
{
    x ^= x >> 12, x ^= x << 25, x ^= x >> 27;
    return x * 0x2545F4914F6CDD1DULL;
}


// Returns the line after the `i`'th one on the pointer chase over `line_c`
// lines, a power of two. It is a full-period linear congruential step, so that
// the chase is a single cycle through every line, in no order that a prefetcher
// could follow.
constexpr uint64_t chase_next(const uint64_t i, const uint64_t line_c)
{
    return (6364136223846793005ULL * i + 1442695040888963407ULL) & (line_c - 1);
}


// Returns the line `k` steps after the `i`'th one on the pointer chase over
// `line_c` lines, composing the chase's step `k` times by squaring.
constexpr uint64_t chase_after(uint64_t i, uint64_t k, const uint64_t line_c)
{
    uint64_t a = 6364136223846793005ULL;    // The composed step is `x → a·x + c`.
    uint64_t c = 1442695040888963407ULL;
    for(; k > 0; k >>= 1)
    {
        if(k & 1)
            i = (a * i + c) & (line_c - 1);

        c = a * c + c, a *= a;
    }

    return i;
}


// Runs `f(t_id, deadline)` with `thread_c` threads until the deadline, which is
// `seconds` away, where `f` returns the number of cache lines that it accessed.
// Counts the events of the group `perf` over the run, if provided, and reports
// the run as the phase `name`. Returns the number of lines accessed and the
// elapsed time in seconds.
template <typename F_>
std::pair<std::size_t, double> run_phase(const std::string& name, const std::size_t thread_c, const double seconds, const F_& f, Perf_Group* const perf)
{
    std::vector<Count> L(thread_c);
    std::vector<std::thread> W;
    W.reserve(thread_c);

    if(perf)
        perf->start();

    const auto t_s = std::chrono::steady_clock::now();
    const auto deadline = t_s + std::chrono::nanoseconds(static_cast<int64_t>(seconds * 1e9));

    for(std::size_t t = 0; t < thread_c; ++t)
        W.emplace_back([&](const std::size_t t_id){ L[t_id].val = f(t_id, deadline); }, t);

    std::for_each(W.begin(), W.end(), [](auto& w){ w.join(); });

    const auto t_e = std::chrono::steady_clock::now();

    if(perf)
        perf->stop();

    std::size_t line_c = 0;
    std::for_each(L.cbegin(), L.cend(), [&](const auto& l){ line_c += l.val; });
    const double t = std::chrono::duration<double>(t_e - t_s).count();
    report_phase(name, line_c, t);
    if(perf)
        report_perf(name, *perf);

    return {line_c, t};
}
