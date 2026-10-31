
#include "bench.hpp"

#include "jellyfish/mer_dna.hpp"
#include "jellyfish/hash_counter.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <thread>
#include <algorithm>
#include <type_traits>


using namespace kache_hash;


// Adapter of Jellyfish's concurrent table, `cooperative::hash_counter`, that
// maps canonical k-mers to counts—a set if `T_` is `void`, inserting its k-mers
// with counts of 0, as PanGenie does. It is configured as `jellyfish count`
// configures it: 7-bit counts that carry into further entries, and up to 126
// reprobes. It doubles in size when full, with all the user threads taking
// part. NB: Jellyfish's table without count bits, as in its bindings' set,
// writes past its end when it doubles, as of v2.3.1.
template <typename T_>
class Jellyfish
{
private:

    typedef jellyfish::cooperative::hash_counter<jellyfish::mer_dna> ht_t;

    static constexpr uint16_t val_len = 7;  // Bits of a count.
    static constexpr uint16_t reprobe_limit = 126;  // Most reprobes of a key.

    // A user thread's k-mers: the current one, its reverse complement, and
    // scratch space for lookups. Each allocates its words, so they are reused.
    struct Mers
    {
        jellyfish::mer_dna fw;
        jellyfish::mer_dna rc;
        jellyfish::mer_dna tmp;
    };

    ht_t M; // The table.
    mutable std::vector<Padded<Mers>> mers; // `mers[t]` are the `t`'th user thread's k-mers.

    // Applies `f` to the canonical form of each of the `n` k-mers of the
    // sequence `s`, rolling them through the k-mers `m`, as Jellyfish's
    // `mer_iterator` does.
    template <typename F_>
    static void for_each_mer(Mers& m, const char* s, std::size_t n, const F_& f);

public:

    static constexpr auto name = "jellyfish";   // Name of the structure.
    static constexpr bool is_set = std::is_void_v<T_>;  // Whether a set.
    static constexpr bool concurrent = true;    // Whether safe for multiple user threads.
    static constexpr bool is_static = false;    // Whether loaded from an index file rather than built.

    // Prints the parameters of the table.
    static void report_params() { print("jellyfish.val_len", val_len); print("jellyfish.reprobe_limit", reprobe_limit); }

    // Constructs a table sized for `init_cap` k-mers, for `thread_c` user
    // threads: with `init_cap` entries, rounded up to a power of 2 by Jellyfish.
    Jellyfish(std::size_t init_cap, std::size_t thread_c);

    // Adds the `n` k-mers of the sequence `s`, for the `t_id`'th user thread.
    void insert(std::size_t t_id, const char* s, std::size_t n);

    // Ends the insertions of the `t_id`'th user thread. Each user thread
    // takes part in doubling the table until all are done.
    void done(std::size_t) { M.done(); }

    // Returns the number of the `n` k-mers of the sequence `s` present in the
    // table, for the `t_id`'th user thread.
    std::size_t find(std::size_t t_id, const char* s, std::size_t n) const;

    // Returns the number of k-mers in the table, counting them in a slice of
    // its entries per user thread.
    std::size_t size() const;

    // Returns the number of entries in the table.
    std::size_t capacity() const { return M.size(); }

    // Prints the statistics of the table.
    void report() const { print("jellyfish.final_reprobe_limit", M.reprobe_limit()); }
};


template <typename T_>
inline Jellyfish<T_>::Jellyfish(const std::size_t init_cap, const std::size_t thread_c):
      M(init_cap, 2 * k, val_len, thread_c, reprobe_limit)
    , mers(thread_c)
{}


template <typename T_>
template <typename F_>
inline void Jellyfish<T_>::for_each_mer(Mers& m, const char* const s, const std::size_t n, const F_& f)
{
    std::size_t filled = 0; // Number of bases rolled into the current k-mer since a non-ACGT base, capped at `k`.
    for(std::size_t i = 0; i < n + k - 1; ++i)
    {
        const int code = m.fw.code(s[i]);
        if(code < 0)
        {
            filled = 0;
            continue;
        }

        m.fw.shift_left(code);
        m.rc.shift_right(m.rc.complement(code));
        filled = std::min<std::size_t>(filled + 1, k);
        if(filled == k)
            f(m.fw < m.rc ? m.fw : m.rc);
    }
}


template <typename T_>
inline void Jellyfish<T_>::insert(const std::size_t t_id, const char* const s, const std::size_t n)
{
    for_each_mer(mers[t_id].unwrap(), s, n, [&](const jellyfish::mer_dna& key)
    {
        if constexpr(is_set)
            M.set(key);
        else
            M.add(key, 1);
    });
}


template <typename T_>
inline std::size_t Jellyfish<T_>::find(const std::size_t t_id, const char* const s, const std::size_t n) const
{
    auto& m = mers[t_id].unwrap();
    const auto& A = *M.ary();   // Doubling replaces the array; it is the last one after all insertions.
    std::size_t h = 0;
    std::size_t id;
    for_each_mer(m, s, n, [&](const jellyfish::mer_dna& key){ h += A.get_key_id(key, &id, m.tmp); });
    return h;
}


template <typename T_>
inline std::size_t Jellyfish<T_>::size() const
{
    const auto& A = *M.ary();
    const auto slice_c = mers.size();
    std::vector<Padded<std::size_t>> C(slice_c, Padded<std::size_t>(0));    // `C[t]` is the number of k-mers in the `t`'th slice.
    std::vector<std::thread> W;
    W.reserve(slice_c);
    for(std::size_t t = 0; t < slice_c; ++t)
        W.emplace_back([&](const std::size_t t_id)
        {
            // Entries holding the high bits of large counts are skipped.
            auto it = A.template iterator_slice<typename ht_t::lazy_iterator>(t_id, slice_c);
            while(it.next())
                C[t_id].unwrap()++;
        }, t);

    std::for_each(W.begin(), W.end(), [](auto& w){ w.join(); });

    std::size_t kmer_c = 0;
    std::for_each(C.cbegin(), C.cend(), [&](const auto& c){ kmer_c += c.unwrap(); });
    return kmer_c;
}


int main(int argc, char* argv[])
{
    const auto cfg = parse_config(argc, argv);

    if(cfg.seq)
        reject("Jellyfish has no sequential variant");

    jellyfish::mer_dna::k(k);   // Jellyfish's k-mers take their length from here.

    if(cfg.set)
        bench<Jellyfish<void>>(cfg);
    else
        bench<Jellyfish<uint64_t>>(cfg);

    return 0;
}
