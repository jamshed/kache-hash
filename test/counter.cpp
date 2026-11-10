
// Checks of the k-mer counting library against exact references. The reader gives the sequences that a sequential
// parser does, for FASTA and FASTQ, of one and of several lines, with "\n" and "\r\n" line ends, from files and pipes,
// in units small enough to split records; and it fails on malformed inputs. The counter's counts, histograms and
// lookups are a brute-force count's. The other mode counts real inputs as PanGenie does, and reports the phases.
// Usage: counter check <dir> [threads]    Writes its inputs into the directory `dir`.
//        counter count <segments> <reads> <threads> <size> [load-factor] [histogram-file]

#include "kache_hash/Fastx_Reader.hpp"
#include "kache_hash/Kmer_Set_Counter.hpp"

#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>


using namespace kache_hash;

constexpr std::size_t k = 31;


// Pseudo-random numbers, from splitmix64.
struct Rng
{
    uint64_t x;

    uint64_t operator()()
    {
        uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }

    std::size_t below(const std::size_t n) { return (*this)() % n; }
};


// A record: its name and its sequence.
typedef std::pair<std::string, std::string> record_t;


// Returns a sequence of length `len` over `alpha`.
static std::string random_seq(Rng& r, const std::size_t len, const std::string& alpha)
{
    std::string s(len, 'A');
    for(auto& c : s)
        c = alpha[r.below(alpha.size())];
    return s;
}


// Returns `c` records of lengths in `[0, max_len]`, but for one in 100 of length `long_len`: mostly bases, in either
// case, with some N's and other IUPAC codes.
static std::vector<record_t> random_records(Rng& r, const std::size_t c, const std::size_t max_len, const std::size_t long_len)
{
    std::vector<record_t> R;
    for(std::size_t i = 0; i < c; ++i)
    {
        const auto len = (r.below(100) == 0 ? long_len : r.below(max_len + 1));
        R.emplace_back("r" + std::to_string(i) + " description", random_seq(r, len, "ACGTACGTACGTACGTACGTacgtNRY"));
    }

    return R;
}


// Returns the FASTA text of the records `R`, with line ends `eol`, but for the last line iff `last_eol`. A sequence is
// in one line if `width` is 0, else in lines of a width of its own in `[1, width]`, with an empty line here and there.
static std::string fasta(const std::vector<record_t>& R, Rng& r, const std::size_t width, const std::string& eol, const bool last_eol = true)
{
    std::string s;
    for(const auto& [name, seq] : R)
    {
        s += ">" + name + eol;
        if(width == 0)
        {
            s += seq + eol;
            continue;
        }

        const auto w = 1 + r.below(width);
        for(std::size_t i = 0; i < seq.size(); i += w)
        {
            s += seq.substr(i, w) + eol;
            if(r.below(50) == 0)
                s += eol;
        }
    }

    if(!last_eol && !s.empty())
        s.resize(s.size() - eol.size());
    return s;
}


// Returns the FASTQ text of the records `R`, as `fasta` does. Some quality lines begin with '@' or '+'.
static std::string fastq(const std::vector<record_t>& R, Rng& r, const std::string& eol, const bool last_eol = true)
{
    std::string s;
    for(const auto& [name, seq] : R)
    {
        auto qual = random_seq(r, seq.size(), "!#%+5:@AEIJ");
        if(!qual.empty() && r.below(4) == 0)
            qual[0] = (r.below(2) ? '@' : '+');
        s += "@" + name + eol + seq + eol + "+" + eol + qual + eol;
    }

    if(!last_eol && !s.empty())
        s.resize(s.size() - eol.size());
    return s;
}


// Returns the non-empty sequences of the FASTA or FASTQ text `s`, as a reader is to give them, sorted, by a sequential
// parse.
static std::vector<std::string> parse(const std::string& s)
{
    std::vector<std::string> L;     // The lines, without their ends.
    for(std::size_t i = 0; i < s.size(); )
    {
        const auto j = std::min(s.find('\n', i), s.size());
        const auto e = (j > i && s[j - 1] == '\r' ? j - 1 : j);
        L.emplace_back(s, i, e - i);
        i = j + 1;
    }

    std::vector<std::string> S;
    if(!s.empty() && s[0] == '>')
    {
        std::string seq;
        for(std::size_t i = 0; i <= L.size(); ++i)
            if(i == L.size() || (!L[i].empty() && L[i][0] == '>'))
            {
                if(!seq.empty())
                    S.push_back(seq);
                seq.clear();
            }
            else
                seq += L[i];
    }
    else if(!s.empty())
        for(std::size_t i = 0; i + 1 < L.size(); i += 4)
            if(!L[i + 1].empty())
                S.push_back(L[i + 1]);

    std::sort(S.begin(), S.end());
    return S;
}


// Writes the text `text` into the file at `path`.
static void write(const std::string& path, const std::string& text)
{
    std::ofstream out(path, std::ios::binary);
    out.write(text.data(), text.size());
}


// Returns the sequences that a reader of `path` gives `thread_c` workers in units of `chunk_bytes` bytes, sorted.
// Clears `ok` if a batch is empty or larger than 1024, or if a worker is not done once.
static std::vector<std::string> read_all(const std::string& path, const std::size_t thread_c, const std::size_t chunk_bytes, bool& ok)
{
    std::vector<std::vector<std::string>> S(thread_c);
    std::vector<int> bad(thread_c, 0), done(thread_c, 0);
    Fastx_Reader(path, chunk_bytes).read(thread_c, [&](const std::size_t t, const char* const* const seqs, const std::size_t* const lens, const std::size_t c)
    {
        bad[t] |= (c == 0 || c > 1024);
        for(std::size_t q = 0; q < c; ++q)
            S[t].emplace_back(seqs[q], lens[q]);
    },
    [&](const std::size_t t){ done[t]++; });

    std::vector<std::string> A;
    for(std::size_t t = 0; t < thread_c; ++t)
    {
        A.insert(A.end(), S[t].begin(), S[t].end());
        ok &= (!bad[t] && done[t] == 1);
    }

    std::sort(A.begin(), A.end());
    return A;
}


// Returns what `read_all` does for the text `text` written by another thread into a pipe at `path`.
static std::vector<std::string> read_pipe(const std::string& path, const std::string& text, const std::size_t thread_c, const std::size_t chunk_bytes, bool& ok)
{
    unlink(path.c_str());
    if(mkfifo(path.c_str(), 0600) != 0)
        throw std::runtime_error("cannot make the pipe " + path);

    std::thread w([&](){ write(path, text); });
    auto S = read_all(path, thread_c, chunk_bytes, ok);
    w.join();
    unlink(path.c_str());
    return S;
}


// Returns whether reading the text `text` from a file at `path`, or from a pipe there iff `pipe`, fails.
static bool fails(const std::string& path, const std::string& text, const bool pipe = false)
{
    if(!pipe)
        write(path, text);
    else
    {
        unlink(path.c_str());
        mkfifo(path.c_str(), 0600);
    }

    std::thread w;
    if(pipe)    // The text fits in the pipe, so that the writer finishes even if the reader stops early.
        w = std::thread([&](){ write(path, text); });

    bool failed = false;
    try
    {
        Fastx_Reader(path, 64).read(3, [](std::size_t, const char* const*, const std::size_t*, std::size_t){});
    }
    catch(const std::runtime_error&)
    {
        failed = true;
    }

    if(pipe)
        w.join();
    return failed;
}


// Checks the reader. Returns the number of failures.
static std::size_t check_reader(const std::string& dir, const std::size_t T)
{
    Rng r{1};
    const auto R = random_records(r, 3000, 400, 70'000);
    const auto Q = random_records(r, 3000, 300, 20'000);
    const std::vector<std::pair<std::string, std::string>> inputs =
    {
        {"fasta",               fasta(R, r, 0, "\n")},
        {"fasta-lines",         fasta(R, r, 80, "\n")},
        {"fasta-lines-crlf",    fasta(R, r, 80, "\r\n")},
        {"fasta-no-last-eol",   fasta(R, r, 0, "\n", false)},
        {"fastq",               fastq(Q, r, "\n")},
        {"fastq-crlf",          fastq(Q, r, "\r\n")},
        {"fastq-no-last-eol",   fastq(Q, r, "\n", false)},
        {"headers",             fasta(random_records(r, 100, 0, 0), r, 0, "\n")},
        {"empty",               ""},
    };

    std::size_t fail = 0;
    for(const auto& [name, text] : inputs)
    {
        const auto path = dir + "/" + name;
        write(path, text);
        const auto exp = parse(text);
        std::size_t runs = 0, bad = 0;
        for(const std::size_t chunk : {std::size_t(1), std::size_t(7), std::size_t(100), std::size_t(4096), std::size_t(1) << 20, Fastx_Reader::chunk_bytes_default})
            for(const auto t : {std::size_t(1), T})
            {
                bool ok = true;
                bad += (read_all(path, t, chunk, ok) != exp || !ok), runs++;
            }

        for(const std::size_t chunk : {std::size_t(16), std::size_t(5000), Fastx_Reader::chunk_bytes_default})
        {
            bool ok = true;
            bad += (read_pipe(path + ".pipe", text, T, chunk, ok) != exp || !ok), runs++;
        }

        std::printf("reader %-18s %8zu bytes, %5zu sequences: %zu runs, %zu wrong\n", name.c_str(), text.size(), exp.size(), runs, bad);
        fail += (bad > 0);
    }

    const std::vector<std::pair<std::string, std::string>> malformed =
    {
        {"short quality",   "@r\nACGT\n+\nII\n"},
        {"no '+' line",     "@r\nACGT\nIIII\n@s\nAC\n+\nII\n"},
        {"cut record",      "@r\nACGT\n+\nIIII\n@s\nAC\n"},
        {"unknown format",  "ACGT\n"},
        {"gzip",            "\x1f\x8b\x08\x00"},
    };

    std::size_t caught = 0;
    for(const auto& [what, text] : malformed)
    {
        const std::size_t c = fails(dir + "/malformed", text) + fails(dir + "/malformed.pipe", text, true);
        if(c < 2)
            std::printf("reader malformed input, %s: read without failing\n", what.c_str());
        caught += c;
    }

    try { Fastx_Reader(dir + "/absent"); } catch(const std::runtime_error&) { caught++; }

    const auto expected = 2 * malformed.size() + 1;
    std::printf("reader malformed inputs: %zu of %zu fail\n", caught, expected);
    fail += (caught != expected);
    return fail;
}


// Returns the canonical 2-bit encoding of the k-mer at `s`, or -1 if it has a non-ACGT character.
static int64_t canonical(const char* const s)
{
    uint64_t f = 0, rc = 0;
    for(std::size_t i = 0; i < k; ++i)
    {
        uint64_t c;
        switch(s[i])
        {
        case 'A': case 'a': c = 0; break;
        case 'C': case 'c': c = 1; break;
        case 'G': case 'g': c = 2; break;
        case 'T': case 't': c = 3; break;
        default: return -1;
        }

        f = (f << 2) | c, rc |= (3 - c) << (2 * i);
    }

    return std::min(f, rc);
}


// Returns the k-mer of the 2-bit encoding `x`.
static std::string decode(const uint64_t x)
{
    std::string s(k, 'A');
    for(std::size_t i = 0; i < k; ++i)
        s[i] = "ACGT"[(x >> (2 * (k - 1 - i))) & 3];
    return s;
}


// Returns the reverse complement of `s`.
static std::string reverse_complement(const std::string& s)
{
    std::string rc(s.rbegin(), s.rend());
    for(auto& c : rc)
        switch(c)
        {
        case 'A': c = 'T'; break; case 'C': c = 'G'; break; case 'G': c = 'C'; break; case 'T': c = 'A'; break;
        case 'a': c = 't'; break; case 'c': c = 'g'; break; case 'g': c = 'c'; break; case 't': c = 'a'; break;
        default: break;
        }

    return rc;
}


// Calls `f(x)` for the canonical encoding `x` of each k-mer of `s` without non-ACGT characters.
template <typename F>
static void for_each_kmer(const std::string& s, const F& f)
{
    for(std::size_t i = 0; i + k <= s.size(); ++i)
        if(const auto x = canonical(s.data() + i); x >= 0)
            f(static_cast<uint64_t>(x));
}


// Checks the counter. Returns the number of failures.
static std::size_t check_counter(const std::string& dir, const std::size_t T)
{
    // The set: the k-mers of a genome, in either case, and of variants of short stretches of it.
    Rng r{2};
    std::vector<record_t> G;
    for(std::size_t i = 0; i < 10; ++i)
    {
        auto g = random_seq(r, 30'000, "ACGTACGTACGTacgt");
        for(std::size_t j = 0; j < 20; ++j)
            g[r.below(g.size())] = 'N';
        G.emplace_back("chr" + std::to_string(i), g);
    }
    for(std::size_t i = 0; i < 2000; ++i)
    {
        const auto& g = G[r.below(10)].second;
        const auto len = 40 + r.below(61);
        auto a = g.substr(r.below(g.size() - len), len);
        a[r.below(len)] = "ACGT"[r.below(4)];
        G.emplace_back("allele" + std::to_string(i), a);
    }

    // The reads: from the set's sequences, on either strand, with errors and N's here and there, or random.
    std::vector<record_t> Q;
    for(std::size_t i = 0; i < 30'000; ++i)
    {
        std::string s;
        if(r.below(4) == 0)
            s = random_seq(r, 100, "ACGT");
        else
        {
            const auto& g = G[r.below(G.size())].second;
            const auto len = std::min<std::size_t>(100, g.size());
            s = g.substr(r.below(g.size() - len + 1), len);
            if(r.below(2))
                s = reverse_complement(s);
            for(std::size_t e = r.below(3); e > 0; --e)
                s[r.below(s.size())] = "ACGTN"[r.below(5)];
        }

        Q.emplace_back("read" + std::to_string(i), s);
    }

    const auto segments = dir + "/segments.fa", reads = dir + "/reads.fq";
    write(segments, fasta(G, r, 60, "\n"));
    const auto reads_text = fastq(Q, r, "\n");
    write(reads, reads_text);

    // The brute-force count.
    std::unordered_map<uint64_t, uint64_t> C;
    for(const auto& g : G)
        for_each_kmer(g.second, [&](const uint64_t x){ C[x] = 0; });
    std::size_t found_exp = 0;
    for(const auto& q : Q)
        for_each_kmer(q.second, [&](const uint64_t x){ if(const auto it = C.find(x); it != C.end()) it->second++, found_exp++; });

    const auto bins = [&](const std::size_t max)
    {
        std::vector<std::size_t> H(max + 1, 0);
        for(const auto& [x, c] : C)
            if(c > 0 && c <= max)
                H[c]++;
        return H;
    };

    std::size_t fail = 0;
    struct Config { std::size_t n, t; bool pipe; };
    for(const auto& [n, t, pipe] : {Config{1000, 1, false}, Config{C.size(), T, false}, Config{C.size(), T, true}})
    {
        Kmer_Set_Counter K(k, n, t);
        K.prime(segments);
        const auto size_primed = K.size();

        std::size_t found;
        if(!pipe)
            found = K.count(reads);
        else
        {
            unlink((reads + ".pipe").c_str()), mkfifo((reads + ".pipe").c_str(), 0600);
            std::thread w([&](){ write(reads + ".pipe", reads_text); });
            found = K.count(reads + ".pipe");
            w.join();
        }

        K.prime(segments);  // Again: the counts stay.

        // Lookups of the set's k-mers, in either orientation and case.
        std::size_t wrong = 0;
        for(const auto& [x, c] : C)
        {
            auto s = decode(x);
            wrong += (K.count_of(s.data()) != c);
            wrong += (K.count_of(reverse_complement(s).data()) != c);
            std::transform(s.begin(), s.end(), s.begin(), [](const char b){ return b | 0x20; });
            wrong += (K.count_of(s.data()) != c);
            s[x % k] = 'N';
            wrong += (K.count_of(s.data()) != 0);
        }

        // Lookups of k-mers absent from the set.
        for(std::size_t i = 0; i < 10'000; ++i)
        {
            const auto s = random_seq(r, k, "ACGT");
            wrong += (!C.count(canonical(s.data())) && K.count_of(s.data()) != 0);
        }

        const bool hist_ok = (K.histogram(3) == bins(3) && K.histogram(10'000) == bins(10'000));
        const bool ok = (size_primed == C.size() && K.size() == C.size() && found == found_exp && !wrong && hist_ok);
        const auto& tm = K.timings();
        std::printf("counter n %6zu, %2zu threads%s: size %zu (ref %zu), found %zu (ref %zu), %zu wrong lookups, histograms %s; "
                    "%.3f s prime, %.3f s count  %s\n",
                    n, t, pipe ? ", reads from a pipe" : "", K.size(), C.size(), found, found_exp, wrong, hist_ok ? "equal" : "DIFFER",
                    tm.prime, tm.count, ok ? "ok" : "<-- WRONG");
        fail += !ok;
    }

    std::size_t caught = 0;
    try { Kmer_Set_Counter(25, 1000, 1); } catch(const std::invalid_argument&) { caught++; }
    try { Kmer_Set_Counter(k, 1000, 0); } catch(const std::invalid_argument&) { caught++; }
    try { Kmer_Set_Counter(k, 1000, 129); } catch(const std::invalid_argument&) { caught++; }
    const bool sup = (Kmer_Set_Counter::supports(31) && !Kmer_Set_Counter::supports(25));
    std::printf("counter unsupported arguments: %zu of 3 fail, supports(31) and not (25): %s\n", caught, sup ? "yes" : "NO");
    fail += (caught != 3 || !sup);
    return fail;
}


// Returns the peak resident set size of the process, in GB.
static double peak_rss_gb()
{
    rusage u;
    getrusage(RUSAGE_SELF, &u);
#ifdef __APPLE__
    return u.ru_maxrss / 1e9;   // Bytes.
#else
    return u.ru_maxrss / 1e6;   // KB.
#endif
}


// Counts the reads' k-mers in the segments' set, as PanGenie does, and reports the phases.
static int count_inputs(const std::string& segments, const std::string& reads, const std::size_t T, const std::size_t n, const double lf, const std::string& hist_path)
{
    Kmer_Set_Counter K(k, n, T, lf);
    K.prime(segments);
    const auto size = K.size();
    const auto found = K.count(reads);
    const auto H = K.histogram(10'000);     // PanGenie's cap.

    std::size_t binned = 0;
    uint64_t digest = 0;
    for(std::size_t c = 0; c < H.size(); ++c)
        binned += H[c], digest = (digest ^ H[c]) * 0x100000001b3ULL;
    if(!hist_path.empty())
    {
        std::ofstream out(hist_path);
        for(std::size_t c = 0; c < H.size(); ++c)
            out << c << '\t' << H[c] << '\n';
    }

    const auto& tm = K.timings();
    std::printf("size=%zu size_after=%zu found=%zu binned=%zu hist_digest=%016llx alloc_s=%.2f prime_s=%.2f count_s=%.2f histogram_s=%.2f peak_rss_gb=%.2f\n",
                size, K.size(), found, binned, static_cast<unsigned long long>(digest), tm.alloc, tm.prime, tm.count, tm.histogram, peak_rss_gb());
    return 0;
}


int main(const int argc, char* argv[])
{
    const std::string mode = (argc > 1 ? argv[1] : "");
    try
    {
        if(mode == "check" && argc >= 3)
        {
            const std::size_t T = (argc > 3 ? std::stoul(argv[3]) : std::max(2u, std::thread::hardware_concurrency()));
            mkdir(argv[2], 0700);
            const auto fail = check_reader(argv[2], T) + check_counter(argv[2], T);
            std::printf("counter checks: %zu failures\n", fail);
            return fail > 0;
        }

        if(mode == "count" && argc >= 6)
            return count_inputs(argv[2], argv[3], std::stoul(argv[4]), static_cast<std::size_t>(std::stod(argv[5])),
                         argc > 6 ? std::stod(argv[6]) : 0.9, argc > 7 ? argv[7] : "");
    }
    catch(const std::exception& e)
    {
        std::fprintf(stderr, "counter: %s\n", e.what());
        return 1;
    }

    std::fprintf(stderr, "Usage: counter check <dir> [threads]\n"
                         "       counter count <segments> <reads> <threads> <size> [load-factor] [histogram-file]\n");
    return 2;
}
