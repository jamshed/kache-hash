
#include "bench.hpp"
#include "kache_hash/Kmer.hpp"

#include <cstdint>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <sys/sysinfo.h>

#include "icebergHT/iceberg_table.h"    // Last: it defines `__restrict__` away for C++.


using namespace kache_hash;


// Adapter of the iceberg table, as a set of k-mers: it can not update the value
// of a present key, and hence can not count.
class Iceberg
{
private:

    static constexpr std::size_t block_slots = (1 << SLOT_BITS) + C_LV2 + MAX_LG_LG_N / D_CHOICES;  // Slots per block, in levels 1 and 2.
    static constexpr double max_lf = 0.96;  // Load factor at which the table resizes, per `RESIZE_THRESHOLD`.

    mutable iceberg_table M;    // The table; reading its size syncs its per-thread counters.
    uint64_t log_slots; // Logarithm of the initial number of level-1 slots.

public:

    static constexpr auto name = "iceberg"; // Name of the structure.
    static constexpr bool is_set = true;    // Whether a set.
    static constexpr bool concurrent = true;    // Whether safe for multiple user threads.
    static constexpr bool is_static = false;    // Whether loaded from an index file rather than built.

    // Prints the parameters of the table.
    static void report_params() {}

    // Constructs a table sized for `init_cap` k-mers, for `thread_c` user
    // threads: the smallest whose levels 1 and 2 hold them below the resize
    // threshold, with at least 8 blocks, as resizes move blocks 8 at a time.
    Iceberg(std::size_t init_cap, std::size_t thread_c);

    // Adds the `n` k-mers of the sequence `s`, for the `t_id`'th user thread.
    void insert(std::size_t t_id, const char* s, std::size_t n);

    // Returns the number of the `n` k-mers of the sequence `s` present in the
    // table.
    std::size_t find(std::size_t t_id, const char* s, std::size_t n) const;

    // Returns the number of k-mers in the table.
    std::size_t size() const { return tot_balls(&M); }

    // Returns the number of k-mer slots in levels 1 and 2 of the table.
    std::size_t capacity() const { return M.metadata.nblocks * block_slots; }

    // Prints the statistics of the table.
    void report() const;
};


inline Iceberg::Iceberg(const std::size_t init_cap, const std::size_t):
    log_slots(std::max<uint64_t>(std::ceil(std::log2(init_cap / (max_lf * block_slots / (1 << SLOT_BITS)))), SLOT_BITS + 3))
{
    if(iceberg_init(&M, log_slots) != 0)
    {
        std::cerr << "Error initializing the iceberg table. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


inline void Iceberg::insert(const std::size_t t_id, const char* const s, const std::size_t n)
{
    for_each_canonical(s, n, [&](const Kmer<k>& key){ iceberg_insert(&M, key.as_int(), 1, static_cast<uint8_t>(t_id)); });
}


inline std::size_t Iceberg::find(const std::size_t t_id, const char* const s, const std::size_t n) const
{
    std::size_t h = 0;
    ValueType val;
    for_each_canonical(s, n, [&](const Kmer<k>& key){ h += iceberg_get_value(&M, key.as_int(), &val, static_cast<uint8_t>(t_id)); });
    return h;
}


inline void Iceberg::report() const
{
    print("iceberg.init_log_slots", log_slots);
    print("iceberg.lv1_size", lv1_balls(&M));
    print("iceberg.lv2_size", lv2_balls(&M));
    print("iceberg.lv3_size", lv3_balls(&M));
    print("table.load_factor", iceberg_load_factor(&M));
}


int main(int argc, char* argv[])
{
    const auto cfg = parse_config(argc, argv);

    if(!cfg.set)
        reject("iceberg can not update a present key's value, and hence can not count");

    if(cfg.seq)
        reject("iceberg has no sequential variant");

    // Its thread ids index its per-thread counters, which are as many as the
    // online processors, unchecked.
    if(cfg.thread_c > static_cast<std::size_t>(get_nprocs()))
        reject("more threads than online processors");

    bench<Iceberg>(cfg);

    return 0;
}
