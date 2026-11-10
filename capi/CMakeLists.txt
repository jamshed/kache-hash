
# The C interface library: concurrent sets of k-mers, and maps from k-mers to values of 8 to 64 bits, for callers in
# other languages, such as Rust. Its interface is plain C, so that projects linking it compile none of the table's code
# and keep their own flags. Each type of table compiles in a source of its own, in parallel.

kache_hash_add_library(kache_capi kache_capi.cpp kache_capi_set.cpp kache_capi_8.cpp kache_capi_16.cpp kache_capi_32.cpp
                       kache_capi_64.cpp)
