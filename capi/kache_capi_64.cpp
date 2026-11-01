
#include "Kmer_Table.hpp"


// The C interface's maps to 64-bit values, of 31-mers with minimizers of length 22.
kache_table* kache_hash::capi::new_map_64(const kache_options& o)
{
    return new Kmer_Table<31, 22, uint64_t>(o);
}
