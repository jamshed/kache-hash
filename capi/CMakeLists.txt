
# The C interface library: concurrent sets of k-mers, and maps from k-mers to
# values of 8 to 64 bits, for callers in other languages, such as Rust. Its
# interface is plain C, so that projects linking it compile none of the
# table's code and keep their own flags. Its own sources take kache-hash's: the
# warning, optimization and ISA flags, and Release's in a build without a type,
# as a parent project's may be. Each type of table compiles in a source of its
# own, in parallel. The library is position-independent, as Rust links
# position-independent executables.

add_library(kache_capi STATIC kache_capi.cpp kache_capi_set.cpp kache_capi_8.cpp kache_capi_16.cpp kache_capi_32.cpp kache_capi_64.cpp)
add_library(kache_hash::kache_capi ALIAS kache_capi)

target_include_directories(kache_capi PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(kache_capi PRIVATE ${COMPILE_FLAGS} $<$<STREQUAL:$<CONFIG>,>:-O3>)
target_compile_definitions(kache_capi PRIVATE $<$<STREQUAL:$<CONFIG>,>:NDEBUG>)
target_link_libraries(kache_capi PRIVATE kache_hash PUBLIC Threads::Threads)
set_target_properties(kache_capi PROPERTIES POSITION_INDEPENDENT_CODE ON)
