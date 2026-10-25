
#include "bench.hpp"
#include "Kmer.hpp"

#include "boost/unordered/concurrent_flat_set.hpp"
#include "boost/unordered/concurrent_flat_map.hpp"
#include "boost/unordered/unordered_flat_set.hpp"
#include "boost/unordered/unordered_flat_map.hpp"

#include <cstdint>
#include <cstddef>
#include <type_traits>


using namespace kache_hash;


// Adapter of the Boost open-addressing tables, that are concurrent iff `mt_`,
// and that map k-mers to `T_` counts—sets if `T_` is `void`.
template <bool mt_, typename T_>
class Boost
{
private:

    typedef std::conditional_t<std::is_void_v<T_>,
        std::conditional_t<mt_, boost::concurrent_flat_set<Kmer<k>, Kmer_Hasher<k>>, boost::unordered_flat_set<Kmer<k>, Kmer_Hasher<k>>>,
        std::conditional_t<mt_, boost::concurrent_flat_map<Kmer<k>, uint64_t, Kmer_Hasher<k>>, boost::unordered_flat_map<Kmer<k>, uint64_t, Kmer_Hasher<k>>>> ht_t;

    ht_t M; // The table.

public:

    static constexpr auto name = "boost";   // Name of the structure.
    static constexpr bool is_set = std::is_void_v<T_>;  // Whether a set.
    static constexpr bool concurrent = mt_; // Whether safe for multiple user threads.
    static constexpr bool is_static = false;    // Whether loaded from an index file rather than built.

    // Prints the parameters of the table.
    static void report_params() { print("boost.container", std::string(mt_ ? "concurrent" : "unordered") + "_flat_" + (is_set ? "set" : "map")); }

    // Constructs a table sized for `init_cap` k-mers, for `thread_c` user
    // threads.
    Boost(std::size_t init_cap, std::size_t thread_c);

    // Adds the `n` k-mers of the sequence `s`, for the `t_id`'th user thread.
    void insert(std::size_t t_id, const char* s, std::size_t n);

    // Returns the number of the `n` k-mers of the sequence `s` present in the
    // table.
    std::size_t find(std::size_t t_id, const char* s, std::size_t n) const;

    // Returns the number of k-mers in the table.
    std::size_t size() const { return M.size(); }

    // Returns the number of k-mer slots in the table.
    std::size_t capacity() const { return M.bucket_count(); }

    // Prints the statistics of the table.
    void report() const { print("table.load_factor", M.load_factor()); }
};


template <bool mt_, typename T_>
inline Boost<mt_, T_>::Boost(const std::size_t init_cap, const std::size_t)
{
    M.reserve(init_cap);
}


template <bool mt_, typename T_>
inline void Boost<mt_, T_>::insert(const std::size_t, const char* const s, const std::size_t n)
{
    for_each_canonical(s, n, [&](const Kmer<k>& key)
    {
        if constexpr(is_set)
            M.insert(key);
        else if constexpr(mt_)
            M.emplace_or_visit(key, 1, [](auto& kv){ kv.second++; });
        else
            M[key]++;
    });
}


template <bool mt_, typename T_>
inline std::size_t Boost<mt_, T_>::find(const std::size_t, const char* const s, const std::size_t n) const
{
    std::size_t h = 0;
    for_each_canonical(s, n, [&](const Kmer<k>& key){ h += M.contains(key); });
    return h;
}


int main(int argc, char* argv[])
{
    const auto cfg = parse_config(argc, argv);

    if(cfg.set)
        cfg.seq ? bench<Boost<false, void>>(cfg) : bench<Boost<true, void>>(cfg);
    else
        cfg.seq ? bench<Boost<false, uint64_t>>(cfg) : bench<Boost<true, uint64_t>>(cfg);

    return 0;
}
