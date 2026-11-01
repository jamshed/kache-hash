
#include "Kmer_Table.hpp"


// The C interface's maps to 8-bit values, of 31-mers with minimizers of length 22.
kache_table* kache_hash::capi::new_map_8(const kache_options& o)
{
    return new Kmer_Table<31, 22, uint8_t>(o);
}
