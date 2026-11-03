
#include "Kmer_Table.hpp"


// The C interface's sets, of 31-mers with minimizers of length 22.
kache_table* kache_hash::capi::new_set(const kache_options& o)
{
    return new Kmer_Table<31, 22, void>(o);
}


bool kache_hash::capi::size_fits(const kache_options& o)
{
    return Streaming_Kmer_Hash_Table<31, true, void, 22>::fits(o.expected_keys, o.load_factor, o.pow2_capacity != 0);
}
