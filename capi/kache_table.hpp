
#ifndef KACHE_HASH_KACHE_TABLE_HPP
#define KACHE_HASH_KACHE_TABLE_HPP



#include "kache_hash/kache_capi.h"

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <mutex>


// A table of the C interface: the operations of a set, or of a map to one
// type of values, behind the C functions of the same names, which validate
// their arguments. Each type's implementation, `kache_hash::capi::Kmer_Table`,
// compiles in a unit of its own.
struct kache_table
{
    const uint32_t value_bytes; // Width of the values; 0 in a set.
    std::mutex reg; // Serializes registrations.
    std::atomic<uint32_t> user_c = 0;   // Number of registered users.

    explicit kache_table(const uint32_t value_bytes): value_bytes(value_bytes) {}

    virtual ~kache_table() = default;

    // Registers the user `user`: the table's next.
    virtual void register_user(uint32_t user) = 0;

    virtual uint64_t insert(uint32_t user, const char* const* seqs, const size_t* lens, size_t n, uint32_t op, uint64_t x) = 0;

    virtual uint64_t insert_oriented(uint32_t user, const char* const* seqs, const size_t* lens, size_t n, uint32_t op, uint64_t x_fwd, uint64_t x_rev) = 0;

    virtual uint64_t update(const char* const* seqs, const size_t* lens, size_t n, uint32_t op, uint64_t x) = 0;

    virtual void count(const char* const* seqs, const size_t* lens, size_t n, uint64_t* counts) const = 0;

    virtual void count_masked(const char* const* seqs, const size_t* lens, size_t n, uint64_t mask_fwd, uint64_t mask_rev, uint64_t* counts) const = 0;

    virtual void sum(const char* const* seqs, const size_t* lens, size_t n, uint64_t* sums) const = 0;

    virtual void lookup(const char* const* seqs, const size_t* lens, size_t n, uint8_t* found, void* values, uint8_t* fwd) const = 0;

    // Returns whether the k-mer at `kmer` is present, with its value at `value`, or 0.
    virtual bool get(const char* kmer, uint64_t& value) const = 0;

    virtual uint64_t size() const = 0;

    virtual uint64_t reduce(uint32_t what, uint64_t mask, uint32_t workers) const = 0;

    virtual void for_each(uint32_t shards, uint32_t shard, kache_keys_fn f, void* ctx) const = 0;

    virtual void get_stats(kache_stats& out) const = 0;
};


namespace kache_hash
{

namespace capi
{

// Returns a new table of 31-mers with the options `o`: a set, or a map to
// values of 8, 16, 32 or 64 bits.
kache_table* new_set(const kache_options& o);
kache_table* new_map_8(const kache_options& o);
kache_table* new_map_16(const kache_options& o);
kache_table* new_map_32(const kache_options& o);
kache_table* new_map_64(const kache_options& o);

// Returns whether a table with the options `o` is within the most buckets
// supported: the same for every type of values.
bool size_fits(const kache_options& o);

}

}



#endif
