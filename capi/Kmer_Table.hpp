
#ifndef KACHE_HASH_KMER_TABLE_HPP
#define KACHE_HASH_KMER_TABLE_HPP



#include "kache_table.hpp"
#include "Streaming_Kmer_Hash_Table.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <limits>
#include <thread>
#include <tuple>
#include <type_traits>
#include <vector>


namespace kache_hash
{

namespace capi
{

// Returns `g(f)`, for the function `f(v)` of the operation `op` with the
// operand `x`.
// NB: Each operation instantiates the table's operations on its own, so that
// none pays for another's in its loops.
template <typename V, typename G>
inline auto with_op(const uint32_t op, const V x, const G& g)
{
    constexpr V max = std::numeric_limits<V>::max();
    switch(op)
    {
    case KACHE_KEEP:    return g([](const V v){ return v; });
    case KACHE_ASSIGN:  return g([x](V){ return x; });
    case KACHE_ADD:     return g([x](const V v){ return static_cast<V>(v > max - x ? max : v + x); });
    case KACHE_OR:      return g([x](const V v){ return static_cast<V>(v | x); });
    case KACHE_MIN:     return g([x](const V v){ return std::min(v, x); });
    default:            return g([x](const V v){ return std::max(v, x); });
    }
}


// Returns `g(f)`, for the function `f(v, fwd)` of the operation `op` with the
// operand `x_fwd` for a k-mer occurring as is, i.e. iff `fwd`, and `x_rev`
// otherwise. The operation is one that an absent k-mer's value, 0, takes to
// the operand.
template <typename V, typename G>
inline auto with_oriented_op(const uint32_t op, const V x_fwd, const V x_rev, const G& g)
{
    constexpr V max = std::numeric_limits<V>::max();
    switch(op)
    {
    case KACHE_ASSIGN:  return g([x_fwd, x_rev](V, const bool fwd){ return fwd ? x_fwd : x_rev; });
    case KACHE_ADD:     return g([x_fwd, x_rev](const V v, const bool fwd){ const V x = (fwd ? x_fwd : x_rev); return static_cast<V>(v > max - x ? max : v + x); });
    case KACHE_OR:      return g([x_fwd, x_rev](const V v, const bool fwd){ return static_cast<V>(v | (fwd ? x_fwd : x_rev)); });
    default:            return g([x_fwd, x_rev](const V v, const bool fwd){ return std::max(v, fwd ? x_fwd : x_rev); });
    }
}


// Runs `f(w)` for each `w` in `[0, c)`, on `c` threads.
template <typename F>
inline void run_on(const uint32_t c, const F& f)
{
    std::vector<std::thread> W;
    W.reserve(c);
    try
    {
        for(uint32_t w = 0; w < c; ++w)
            W.emplace_back(f, w);
    }
    catch(...) // The threads started are to be joined before their handles go.
    {
        for(auto& t : W)
            t.join();
        throw;
    }

    for(auto& t : W)
        t.join();
}


// A table of the C interface: of `k`-mers, with minimizers of length `l`, and
// values of type `V`; a set iff `V` is `void`.
template <uint16_t k, uint16_t l, typename V>
class Kmer_Table final : public kache_table
{
private:

    static constexpr bool is_set = std::is_void_v<V>;
    typedef std::conditional_t<is_set, uint8_t, V> val_t;   // The values' type; a placeholder in a set.
    typedef Streaming_Kmer_Hash_Table<k, true, V, l> table_t;

    table_t M;
    std::array<typename table_t::Token, KACHE_MAX_USERS> token;   // The registered users'.

public:

    explicit Kmer_Table(const kache_options& o):
          kache_table(o.value_bytes)
        , M(o.expected_keys, o.resize_workers, o.load_factor, o.pow2_capacity != 0)
    {}

    void register_user(const uint32_t user) override { token[user] = M.register_user(); }

    uint64_t insert(const uint32_t user, const char* const* const seqs, const size_t* const lens, const size_t n, const uint32_t op, const uint64_t x) override
    {
        if constexpr(is_set)
            return (void)op, (void)x, M.insert(seqs, lens, n, token[user]);
        else
            return with_op(op, static_cast<V>(x), [&](const auto& f){ return M.upsert(seqs, lens, n, f, static_cast<V>(x), token[user]); });
    }

    uint64_t insert_oriented(const uint32_t user, const char* const* const seqs, const size_t* const lens, const size_t n, const uint32_t op, const uint64_t x_fwd, const uint64_t x_rev) override
    {
        if constexpr(is_set)
            return (void)user, (void)seqs, (void)lens, (void)n, (void)op, (void)x_fwd, (void)x_rev, 0;
        else
            return with_oriented_op(op, static_cast<V>(x_fwd), static_cast<V>(x_rev), [&](const auto& f){ return M.upsert_oriented(seqs, lens, n, f, token[user]); });
    }

    uint64_t update(const char* const* const seqs, const size_t* const lens, const size_t n, const uint32_t op, const uint64_t x) override
    {
        if constexpr(is_set)
            return (void)seqs, (void)lens, (void)n, (void)op, (void)x, 0;
        else
            return with_op(op, static_cast<V>(x), [&](const auto& f){ return M.update(seqs, lens, n, f); });
    }

    void count(const char* const* const seqs, const size_t* const lens, const size_t n, uint64_t* const counts) const override
    {
        std::fill(counts, counts + n, 0);
        M.find(seqs, lens, n, [&](const std::size_t q, std::size_t, const auto& r){ counts[q] += static_cast<bool>(r); });
    }

    void count_masked(const char* const* const seqs, const size_t* const lens, const size_t n, const uint64_t mask_fwd, const uint64_t mask_rev, uint64_t* const counts) const override
    {
        if constexpr(is_set)
            (void)seqs, (void)lens, (void)n, (void)mask_fwd, (void)mask_rev, (void)counts;
        else
        {
            const auto m_fwd = static_cast<V>(mask_fwd), m_rev = static_cast<V>(mask_rev);
            std::fill(counts, counts + n, 0);
            M.find(seqs, lens, n, [&](const std::size_t q, std::size_t, const auto& r, const bool fwd){ counts[q] += (r && (*r & (fwd ? m_fwd : m_rev)) != 0); });
        }
    }

    void sum(const char* const* const seqs, const size_t* const lens, const size_t n, uint64_t* const sums) const override
    {
        if constexpr(is_set)
            (void)seqs, (void)lens, (void)n, (void)sums;
        else
        {
            std::fill(sums, sums + n, 0);
            M.find(seqs, lens, n, [&](const std::size_t q, std::size_t, const auto& r){ sums[q] += r.value_or(0); });
        }
    }

    void lookup(const char* const* const seqs, const size_t* const lens, const size_t n, uint8_t* const found, void* const values, uint8_t* const fwd) const override
    {
        std::vector<uint64_t> o(n);  // Index of each sequence's first position in the outputs.
        uint64_t p = 0;
        for(std::size_t i = 0; i < n; ++i)
            o[i] = p, p += (lens[i] >= k ? lens[i] - k + 1 : 0);

        // Positions with a non-ACGT character get no callback.
        val_t* const val = static_cast<val_t*>(values);
        std::fill(found, found + p, uint8_t(KACHE_NOT_KMER));
        if(val != nullptr)
            std::fill(val, val + p, 0);
        if(fwd != nullptr)
            std::fill(fwd, fwd + p, 0);

        const auto put = [&](const uint64_t j, const auto& r)
        {
            found[j] = (r ? KACHE_PRESENT : KACHE_ABSENT);
            if constexpr(!is_set)
                if(val != nullptr)
                    val[j] = r.value_or(0);
        };

        if(fwd == nullptr)
            M.find(seqs, lens, n, [&](const std::size_t q, const std::size_t i, const auto& r){ put(o[q] + i, r); });
        else
            M.find(seqs, lens, n, [&](const std::size_t q, const std::size_t i, const auto& r, const bool f){ put(o[q] + i, r), fwd[o[q] + i] = f; });
    }

    bool get(const char* const kmer, uint64_t& value) const override
    {
        const auto r = M.find(kmer);
        if constexpr(is_set)
            return value = 0, r;
        else
            return value = r.value_or(0), r.has_value();
    }

    uint64_t size() const override { return M.size(); }

    uint64_t reduce(const uint32_t what, const uint64_t mask, const uint32_t workers) const override
    {
        std::vector<uint64_t> part(workers, 0);
        run_on(workers, [&](const uint32_t w)
        {
            uint64_t c = 0;
            if constexpr(is_set)
                (void)what, (void)mask, M.for_each([&](const Kmer<k>&){ c++; }, workers, w);
            else
            {
                const auto m = static_cast<V>(mask);
                if(what == KACHE_REDUCE_KEYS)
                    M.for_each([&](const Kmer<k>&, const V v){ c += ((v & m) != 0); }, workers, w);
                else if(what == KACHE_REDUCE_POPCOUNT)
                    M.for_each([&](const Kmer<k>&, const V v){ c += std::popcount(static_cast<V>(v & m)); }, workers, w);
                else
                    M.for_each([&](const Kmer<k>&, const V v){ c += static_cast<V>(v & m); }, workers, w);
            }

            part[w] = c;
        });

        uint64_t r = 0;
        for(const auto c : part)
            r += c;

        return r;
    }

    void for_each(const uint32_t shards, const uint32_t shard, const kache_keys_fn f, void* const ctx) const override
    {
        constexpr std::size_t chunk_sz = 4096;  // Most keys per call of `f`.
        std::vector<uint64_t> K(chunk_sz);
        std::vector<val_t> W(is_set ? 0 : chunk_sz);
        std::size_t c = 0;  // Number of keys in the chunk.
        bool done = false;  // Whether `f` asked for no more.
        const auto flush = [&]()
        {
            if(c > 0 && !done)
                done = (f(ctx, K.data(), is_set ? nullptr : W.data(), c) != 0);
            c = 0;
        };

        if constexpr(is_set)
            M.for_each([&](const Kmer<k>& x)
            {
                if(done)
                    return;

                K[c++] = x.as_int();
                if(c == chunk_sz)
                    flush();
            }, shards, shard);
        else
            M.for_each([&](const Kmer<k>& x, const V v)
            {
                if(done)
                    return;

                K[c] = x.as_int(), W[c] = v;
                if(++c == chunk_sz)
                    flush();
            }, shards, shard);

        flush();
    }

    void get_stats(kache_stats& out) const override
    {
        out = kache_stats{};
        out.capacity = M.capacity();
        out.bytes = M.bucket_bytes();
        out.keys = M.size();
        out.overflow_keys = M.overflow_size();
        if(out.keys > 0)    // The distribution of an empty table is undefined.
            std::tie(out.primary_pct, out.secondary_pct, out.overflow_pct) = M.key_distribution();
    }
};

}

}



#endif
