
#include "kache_capi.h"
#include "kache_table.hpp"

#include <exception>
#include <mutex>
#include <new>


namespace
{

// Runs `f`, returning its status, or that of the exception it throws.
template <typename F>
int guard(const F& f) noexcept
{
    try
    {
        return f();
    }
    catch(const std::bad_alloc&)
    {
        return KACHE_ENOMEM;
    }
    catch(...)
    {
        return KACHE_EINTERNAL;
    }
}


// Returns whether the sequences `seqs[0..n)` of lengths `lens[0..n)` are
// valid: only the empty ones may be NULL.
bool valid(const char* const* const seqs, const size_t* const lens, const size_t n)
{
    if(n == 0)
        return true;
    if(seqs == nullptr || lens == nullptr)
        return false;

    for(size_t i = 0; i < n; ++i)
        if(seqs[i] == nullptr && lens[i] > 0)
            return false;

    return true;
}


// Returns whether `t` is a map, of values that the operand `x` fits.
bool fits(const kache_table* const t, const uint64_t x)
{
    return t->value_bytes > 0 && (t->value_bytes == 8 || x >> (8 * t->value_bytes) == 0);
}


// Returns whether `user` is a registered user of the table `t`.
bool registered(const kache_table* const t, const uint32_t user)
{
    return user < t->user_c.load(std::memory_order_acquire);
}

}


extern "C"
{

int kache_supports(const uint32_t k, const uint32_t value_bytes) noexcept
{
    return k == 31 && (value_bytes == 0 || value_bytes == 1 || value_bytes == 2 || value_bytes == 4 || value_bytes == 8);
}


kache_options kache_default_options(const uint32_t k, const uint32_t value_bytes) noexcept
{
    return kache_options{k, value_bytes, uint64_t(1) << 22, 0.8, 16, 1};
}


int kache_new(const kache_options* const options, kache_table** const out) noexcept
{
    if(out == nullptr)
        return KACHE_EINVAL;

    *out = nullptr;
    if(options == nullptr)
        return KACHE_EINVAL;

    const auto& o = *options;
    if(!kache_supports(o.k, o.value_bytes) || !(o.load_factor > 0 && o.load_factor < 1) || o.resize_workers < 1 ||
       o.pow2_capacity > 1 || !kache_hash::capi::size_fits(o))
        return KACHE_EINVAL;

    return guard([&]()
    {
        using namespace kache_hash::capi;
        switch(o.value_bytes)
        {
        case 0:     *out = new_set(o);      break;
        case 1:     *out = new_map_8(o);    break;
        case 2:     *out = new_map_16(o);   break;
        case 4:     *out = new_map_32(o);   break;
        default:    *out = new_map_64(o);   break;
        }

        return KACHE_OK;
    });
}


void kache_free(kache_table* const t) noexcept
{
    delete t;
}


int kache_register(kache_table* const t, uint32_t* const user) noexcept
{
    if(t == nullptr || user == nullptr)
        return KACHE_EINVAL;

    return guard([&]()
    {
        // The table doesn't bound its users, so the interface does.
        const std::lock_guard<std::mutex> lock(t->reg);
        const auto u = t->user_c.load(std::memory_order_relaxed);
        if(u == KACHE_MAX_USERS)
            return KACHE_EUSERS;

        t->register_user(u);
        t->user_c.store(u + 1, std::memory_order_release);
        *user = u;
        return KACHE_OK;
    });
}


int kache_insert(kache_table* const t, const uint32_t user, const char* const* const seqs, const size_t* const lens, const size_t n, const uint32_t op, const uint64_t x, uint64_t* const inserted) noexcept
{
    if(t == nullptr || !registered(t, user) || !valid(seqs, lens, n) || (t->value_bytes > 0 && (op > KACHE_MAX || !fits(t, x))))
        return KACHE_EINVAL;

    return guard([&]()
    {
        const auto c = t->insert(user, seqs, lens, n, op, x);
        if(inserted != nullptr)
            *inserted = c;

        return KACHE_OK;
    });
}


int kache_insert_oriented(kache_table* const t, const uint32_t user, const char* const* const seqs, const size_t* const lens, const size_t n, const uint32_t op, const uint64_t x_fwd, const uint64_t x_rev, uint64_t* const inserted) noexcept
{
    const bool op_ok = (op == KACHE_ASSIGN || op == KACHE_ADD || op == KACHE_OR || op == KACHE_MAX);
    if(t == nullptr || !registered(t, user) || !valid(seqs, lens, n) || !op_ok || !fits(t, x_fwd) || !fits(t, x_rev))
        return KACHE_EINVAL;

    return guard([&]()
    {
        const auto c = t->insert_oriented(user, seqs, lens, n, op, x_fwd, x_rev);
        if(inserted != nullptr)
            *inserted = c;

        return KACHE_OK;
    });
}


int kache_update(kache_table* const t, const char* const* const seqs, const size_t* const lens, const size_t n, const uint32_t op, const uint64_t x, uint64_t* const present) noexcept
{
    if(t == nullptr || !valid(seqs, lens, n) || op > KACHE_MAX || !fits(t, x))
        return KACHE_EINVAL;

    return guard([&]()
    {
        const auto c = t->update(seqs, lens, n, op, x);
        if(present != nullptr)
            *present = c;

        return KACHE_OK;
    });
}


int kache_count(const kache_table* const t, const char* const* const seqs, const size_t* const lens, const size_t n, uint64_t* const counts) noexcept
{
    if(t == nullptr || !valid(seqs, lens, n) || (n > 0 && counts == nullptr))
        return KACHE_EINVAL;

    return guard([&](){ return t->count(seqs, lens, n, counts), KACHE_OK; });
}


int kache_count_masked(const kache_table* const t, const char* const* const seqs, const size_t* const lens, const size_t n, const uint64_t mask_fwd, const uint64_t mask_rev, uint64_t* const counts) noexcept
{
    if(t == nullptr || t->value_bytes == 0 || !valid(seqs, lens, n) || (n > 0 && counts == nullptr))
        return KACHE_EINVAL;

    return guard([&](){ return t->count_masked(seqs, lens, n, mask_fwd, mask_rev, counts), KACHE_OK; });
}


int kache_sum(const kache_table* const t, const char* const* const seqs, const size_t* const lens, const size_t n, uint64_t* const sums) noexcept
{
    if(t == nullptr || t->value_bytes == 0 || !valid(seqs, lens, n) || (n > 0 && sums == nullptr))
        return KACHE_EINVAL;

    return guard([&](){ return t->sum(seqs, lens, n, sums), KACHE_OK; });
}


int kache_lookup(const kache_table* const t, const char* const* const seqs, const size_t* const lens, const size_t n, uint8_t* const found, void* const values, uint8_t* const fwd) noexcept
{
    if(t == nullptr || !valid(seqs, lens, n) || (n > 0 && found == nullptr) || (t->value_bytes == 0 && values != nullptr))
        return KACHE_EINVAL;

    return guard([&](){ return t->lookup(seqs, lens, n, found, values, fwd), KACHE_OK; });
}


int kache_get(const kache_table* const t, const char* const kmer, uint8_t* const found, uint64_t* const value) noexcept
{
    if(t == nullptr || kmer == nullptr || found == nullptr)
        return KACHE_EINVAL;

    return guard([&]()
    {
        uint64_t v;
        *found = t->get(kmer, v);
        if(value != nullptr)
            *value = v;

        return KACHE_OK;
    });
}


int kache_size(const kache_table* const t, uint64_t* const keys) noexcept
{
    if(t == nullptr || keys == nullptr)
        return KACHE_EINVAL;

    return guard([&](){ return *keys = t->size(), KACHE_OK; });
}


int kache_reduce(const kache_table* const t, const uint32_t what, const uint64_t mask, const uint32_t workers, uint64_t* const out) noexcept
{
    if(t == nullptr || out == nullptr || workers < 1 || what > KACHE_REDUCE_SUM || (t->value_bytes == 0 && what != KACHE_REDUCE_KEYS))
        return KACHE_EINVAL;

    return guard([&](){ return *out = t->reduce(what, mask, workers), KACHE_OK; });
}


int kache_for_each(const kache_table* const t, const uint32_t shards, const uint32_t shard, const kache_keys_fn f, void* const ctx) noexcept
{
    if(t == nullptr || f == nullptr || shard >= shards)
        return KACHE_EINVAL;

    return guard([&](){ return t->for_each(shards, shard, f, ctx), KACHE_OK; });
}


int kache_get_stats(const kache_table* const t, kache_stats* const out) noexcept
{
    if(t == nullptr || out == nullptr)
        return KACHE_EINVAL;

    return guard([&](){ return t->get_stats(*out), KACHE_OK; });
}


const char* kache_strerror(const int code) noexcept
{
    switch(code)
    {
    case KACHE_OK:          return "success";
    case KACHE_EINVAL:      return "invalid argument";
    case KACHE_EUSERS:      return "too many users registered";
    case KACHE_ENOMEM:      return "out of memory";
    case KACHE_EINTERNAL:   return "internal error";
    default:                return "unknown status";
    }
}

}
