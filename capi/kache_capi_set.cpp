
#include "Kmer_Table.hpp"


// The C interface's sets, of 31-mers with minimizers of length 22.
kache_table* kache_hash::capi::new_set(const kache_options& o)
{
    return new Kmer_Table<31, 22, void>(o);
}
