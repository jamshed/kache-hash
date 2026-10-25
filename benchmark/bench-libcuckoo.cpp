
#include "bench.hpp"
#include "Kmer.hpp"

#include "libcuckoo/cuckoohash_map.hh"

#include <cstdint>
#include <cstddef>
#include <type_traits>


using namespace kache_hash;


// Adapter of the libcuckoo table, that maps k-mers to `T_` counts—a set if `T_`
// is `void`, stored as a map to `bool`s.
template <typename T_>
class Libcuckoo
{
private:

    typedef libcuckoo::cuckoohash_map<Kmer<k>, std::conditional_t<std::is_void_v<T_>, bool, uint64_t>, Kmer_Hasher<k>> ht_t;

    ht_t M; // The table.

public:

    static constexpr auto name = "libcuckoo";   // Name of the structure.
    static constexpr bool is_set = std::is_void_v<T_>;  // Whether a set.
    static constexpr bool concurrent = true;    // Whether safe for multiple user threads.

    // Prints the parameters of the table.
    static void report_params() {}

    // Constructs a table sized for `init_cap` k-mers, for `thread_c` user
    // threads.
    Libcuckoo(std::size_t init_cap, std::size_t thread_c);

    // Adds the `n` k-mers of the sequence `s`, for the `t_id`'th user thread.
    void insert(std::size_t t_id, const char* s, std::size_t n);

    // Returns the number of the `n` k-mers of the sequence `s` present in the
    // table.
    std::size_t find(std::size_t t_id, const char* s, std::size_t n) const;

    // Returns the number of k-mers in the table.
    std::size_t size() const { return M.size(); }

    // Returns the number of k-mer slots in the table.
    std::size_t capacity() const { return M.capacity(); }

    // Prints the statistics of the table.
    void report() const { print("table.load_factor", M.load_factor()); }
};


template <typename T_>
inline Libcuckoo<T_>::Libcuckoo(const std::size_t init_cap, const std::size_t):
    M(init_cap)
{}


template <typename T_>
inline void Libcuckoo<T_>::insert(const std::size_t, const char* const s, const std::size_t n)
{
    for_each_canonical(s, n, [&](const Kmer<k>& key)
    {
        if constexpr(is_set)
            M.insert(key, true);
        else
            M.upsert(key, [](uint64_t& c){ c++; }, 1);
    });
}


template <typename T_>
inline std::size_t Libcuckoo<T_>::find(const std::size_t, const char* const s, const std::size_t n) const
{
    std::size_t h = 0;
    for_each_canonical(s, n, [&](const Kmer<k>& key){ h += M.contains(key); });
    return h;
}


int main(int argc, char* argv[])
{
    const auto cfg = parse_config(argc, argv);

    if(cfg.seq)
        reject("libcuckoo has no sequential variant");

    if(cfg.set)
        bench<Libcuckoo<void>>(cfg);
    else
        bench<Libcuckoo<uint64_t>>(cfg);

    return 0;
}
