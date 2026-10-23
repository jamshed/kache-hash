
#include "bench.hpp"
#include "Streaming_Kmer_Hash_Table.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <type_traits>
#include <algorithm>


using namespace kache_hash;


// Adapter of the kache-hash table, that is concurrent iff `mt_`, and that maps
// k-mers to `T_` counts—a set if `T_` is `void`.
template <bool mt_, typename T_>
class Kache_Hash
{
private:

    typedef Streaming_Kmer_Hash_Table<k, mt_, T_, l> ht_t;

    ht_t M; // The table.
    std::vector<typename ht_t::Token> token;    // `token[t]` is the token of the `t`'th user thread.

public:

    static constexpr auto name = "kache-hash";  // Name of the structure.
    static constexpr bool is_set = std::is_void_v<T_>;  // Whether a set.
    static constexpr bool concurrent = mt_; // Whether safe for multiple user threads.

    // Prints the parameters of the table.
    static void report_params() { print("l", l); }

    // Constructs a table sized for `init_cap` k-mers, for `thread_c` user
    // threads.
    Kache_Hash(std::size_t init_cap, std::size_t thread_c);

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
    void report() const;
};


template <bool mt_, typename T_>
inline Kache_Hash<mt_, T_>::Kache_Hash(const std::size_t init_cap, const std::size_t thread_c):
      M(init_cap, thread_c)
    , token(thread_c)
{
    std::for_each(token.begin(), token.end(), [&](auto& t){ t = M.register_user(); });
}


template <bool mt_, typename T_>
inline void Kache_Hash<mt_, T_>::insert(const std::size_t t_id, const char* const s, const std::size_t n)
{
    for_each_window(s, n, [&](const Kmer_Window<k, l>& w)
    {
        if constexpr(is_set)
            M.insert(w, token[t_id]);
        else
            M.upsert(w, [](const uint64_t v){ return v + 1; }, 1, token[t_id]);
    });
}


template <bool mt_, typename T_>
inline std::size_t Kache_Hash<mt_, T_>::find(const std::size_t, const char* const s, const std::size_t n) const
{
    std::size_t h = 0;
    for_each_window(s, n, [&](const Kmer_Window<k, l>& w){ h += static_cast<bool>(M.find(w)); });
    return h;
}


template <bool mt_, typename T_>
inline void Kache_Hash<mt_, T_>::report() const
{
    print("table.overflow_size", M.overflow_size());

    const auto [primary, secondary, overflow] = M.key_distribution();
    print("table.primary_pct", primary);
    print("table.secondary_pct", secondary);
    print("table.overflow_pct", overflow);
}


int main(int argc, char* argv[])
{
    const auto cfg = parse_config(argc, argv);

    // The table supports up to 64 users.
    if(cfg.thread_c > 64)
        reject("more than 64 threads");

    if(cfg.set)
        cfg.seq ? bench<Kache_Hash<false, void>>(cfg) : bench<Kache_Hash<true, void>>(cfg);
    else
        cfg.seq ? bench<Kache_Hash<false, uint64_t>>(cfg) : bench<Kache_Hash<true, uint64_t>>(cfg);

    return 0;
}
