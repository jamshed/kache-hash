
#ifndef KACHE_HASH_HARNESS_HPP
#define KACHE_HASH_HARNESS_HPP



#include "Perf_Group.hpp"
#include "Parser.hpp"
#include "DNA_Utility.hpp"
#include "utility.hpp"

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <thread>
#include <algorithm>

#include <sys/resource.h>


using namespace kache_hash;


// The minimizer size of the kache-hash tables is 19, unless a build sets `MINIMIZER_SIZE`, as the drivers of the
// minimizer-size sweep do.
#ifndef MINIMIZER_SIZE
#define MINIMIZER_SIZE 19
#endif

constexpr uint16_t k = 31;  // k-mer size.
constexpr uint16_t l = MINIMIZER_SIZE;  // Minimizer size of the kache-hash tables.


// Fragments of the sequences in an input: their maximal stretches of at least
// `k` bases without placeholders.
class Fragments
{
private:

    std::string seq;    // Concatenated fragments.
    std::vector<std::size_t> off;   // `off[i]` is the offset of the `i`'th fragment in `seq`; the last entry is the size of `seq`.
    std::size_t seq_c;  // Number of sequences in the input.
    std::size_t base_c; // Number of bases in the input.

public:

    // Reads the fragments of the input at path `path`.
    Fragments(const std::string& path);

    // Returns the number of fragments.
    std::size_t size() const { return off.size() - 1; }

    // Returns the `i`'th fragment.
    const char* operator[](const std::size_t i) const { return seq.data() + off[i]; }

    // Returns the number of k-mers in the `i`'th fragment.
    std::size_t kmer_count(const std::size_t i) const { return off[i + 1] - off[i] - (k - 1); }

    // Returns the number of k-mers in the fragments.
    std::size_t kmer_count() const { return seq.size() - size() * (k - 1); }

    // Returns the number of sequences in the input.
    std::size_t seq_count() const { return seq_c; }

    // Returns the number of bases in the input.
    std::size_t base_count() const { return base_c; }

    // Returns the number of bases in the input that are not in any fragment.
    std::size_t skipped_base_count() const { return base_c - seq.size(); }
};


inline Fragments::Fragments(const std::string& path):
      seq_c(0)
    , base_c(0)
{
    Parser input(path);
    while(input.read_next_seq())
    {
        const auto s = input.seq();
        const auto len = input.seq_len();
        seq_c++, base_c += len;

        std::size_t start = 0;  // Start of the current stretch without placeholders.
        for(std::size_t i = 0; i <= len; ++i)
            if(i == len || DNA_Utility::is_placeholder(s[i]))
            {
                if(i - start >= k)
                    off.push_back(seq.size()), seq.append(s + start, i - start);

                start = i + 1;
            }
    }

    input.close();
    off.push_back(seq.size());
}


// Prints the measurement `key` with value `val`.
template <typename T_>
void print(const std::string& key, const T_& val)
{
    std::cout << key << "\t" << val << std::endl;
}


// Prints the statistics of the input `name` at path `path`, having the
// fragments `F`.
inline void report_input(const std::string& name, const std::string& path, const Fragments& F)
{
    print(name, path);
    print(name + ".sequences", F.seq_count());
    print(name + ".bases", F.base_count());
    print(name + ".skipped_bases", F.skipped_base_count());
    print(name + ".fragments", F.size());
    print(name + ".kmers", F.kmer_count());
}


// Prints the measurements of the phase `name` that performed `op_c` operations
// in `t` seconds.
inline void report_phase(const std::string& name, const std::size_t op_c, const double t)
{
    print(name + ".ops", op_c);
    print(name + ".seconds", t);
    print(name + ".mops", op_c / (t * 1'000'000.0));
}


// Prints the counts of the performance-counter events in the group `perf` for
// the phase `name`.
inline void report_perf(const std::string& name, const Perf_Group& perf)
{
    for(std::size_t i = 0; i < perf.size(); ++i)
        print("perf." + name + "." + perf.name(i), perf.count(i));

    print("perf." + name + ".running_pct", perf.running_pct());
}


// Returns the maximum resident set size of the process so far, in bytes.
inline std::size_t max_rss_bytes()
{
#ifdef __APPLE__
    constexpr std::size_t rss_unit = 1; // `ru_maxrss` is in bytes on macOS, and in kilobytes on Linux.
#else
    constexpr std::size_t rss_unit = 1024;
#endif

    rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return ru.ru_maxrss * rss_unit;
}


// Prints the maximum resident set size of the process so far, as of the point
// `name` in the run.
inline void report_rss(const std::string& name)
{
    print("memory." + name + "_max_rss_bytes", max_rss_bytes());
}


// Prints the anonymous memory of the process, and how much of it transparent
// huge pages back, as of the point `name` in the run.
inline void report_huge_pages(const std::string& name)
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
#else
    (void)name;
#endif
}


// Prints the maximum resident set size of the process, and the mode of its
// transparent huge pages.
inline void report_memory()
{
    print("memory.max_rss_bytes", max_rss_bytes());

#ifdef __linux__
    std::ifstream thp("/sys/kernel/mm/transparent_hugepage/enabled");   // Lists the modes, e.g. "always [madvise] never".
    std::string mode;
    while(thp >> mode)
        if(mode.front() == '[')
            print("memory.thp", mode.substr(1, mode.size() - 2));
#endif
}


// Runs `f(t_id, s, n)` in parallel with `thread_c` threads, where the `t_id`'th
// thread gets the `t_id`'th contiguous share of the k-mers in the fragments
// `F`, in runs of `n` k-mers starting at `s`, each within a fragment. Counts the
// events of the group `perf` over the run, if provided. Returns the elapsed time
// in seconds.
template <typename F_>
double run_parallel(const Fragments& F, const std::size_t thread_c, const F_& f, Perf_Group* const perf = nullptr)
{
    const auto kmer_c = F.kmer_count();
    const auto share_start = [&](const std::size_t t){ return t * kmer_c / thread_c; };

    std::vector<std::pair<std::size_t, std::size_t>> first(thread_c);   // The `t`'th share starts at the `first[t].second`'th k-mer of the `first[t].first`'th fragment.
    std::size_t frag = 0;
    std::size_t pre = 0;    // Number of k-mers in the fragments before the `frag`'th.
    for(std::size_t t = 0; t < thread_c; ++t)
    {
        while(frag < F.size() && pre + F.kmer_count(frag) <= share_start(t))
            pre += F.kmer_count(frag++);

        first[t] = {frag, share_start(t) - pre};
    }

    std::vector<std::thread> W;
    W.reserve(thread_c);

    if(perf)
        perf->start();

    const auto t_s = now();

    for(std::size_t t = 0; t < thread_c; ++t)
        W.emplace_back([&](const std::size_t t_id)
        {
            auto [i, p] = first[t_id];
            for(auto n = share_start(t_id + 1) - share_start(t_id); n > 0; i++, p = 0)
            {
                const auto m = std::min(n, F.kmer_count(i) - p);
                f(t_id, F[i] + p, m);
                n -= m;
            }
        }, t);

    std::for_each(W.begin(), W.end(), [](auto& w){ w.join(); });

    const auto t_e = now();

    if(perf)
        perf->stop();

    return duration(t_e - t_s);
}


// Runs `f(t_id, seqs, lens, c)` in parallel with `thread_c` threads, where the
// `t_id`'th thread gets the `t_id`'th contiguous share of the k-mers in the
// fragments `F` as in `run_parallel`, in batches of up to `batch_sz` pieces of
// fragments: `seqs[i]` and `lens[i]` are a piece's start and its length in
// bases. Counts the events of the group `perf` over the run, if provided.
// Returns the elapsed time in seconds.
template <typename F_>
double run_parallel_batched(const Fragments& F, const std::size_t thread_c, const std::size_t batch_sz, const F_& f, Perf_Group* const perf = nullptr)
{
    struct Batch
    {
        std::vector<const char*> seqs;
        std::vector<std::size_t> lens;
    };

    std::vector<Padded<Batch>> B(thread_c);   // Padded, as each thread updates its own often.
    const auto flush = [&](const std::size_t t_id)
    {
        auto& b = B[t_id].unwrap();
        if(!b.seqs.empty())
            f(t_id, b.seqs.data(), b.lens.data(), b.seqs.size());
        b.seqs.clear(), b.lens.clear();
    };

    // The pieces of a thread's share come in order; the thread's last piece ends its share.
    const auto kmer_c = F.kmer_count();
    std::vector<Padded<std::size_t>> left(thread_c);    // Number of k-mers left in each share.
    for(std::size_t t = 0; t < thread_c; ++t)
        left[t].unwrap() = (t + 1) * kmer_c / thread_c - t * kmer_c / thread_c;

    return run_parallel(F, thread_c, [&](const std::size_t t_id, const char* const s, const std::size_t n)
    {
        auto& b = B[t_id].unwrap();
        b.seqs.push_back(s), b.lens.push_back(n + k - 1);
        left[t_id].unwrap() -= n;
        if(b.seqs.size() == batch_sz || left[t_id].unwrap() == 0)
            flush(t_id);
    }, perf);
}



#endif
