
# Consumers of kache-hash, for test/consumer/run.sh: the header-only table, the C interface, and the counter, from an
# installed kache-hash that `find_package` finds, or from its source tree through FetchContent, which
# FIND_PACKAGE_ARGS lets take an installed one first.

cmake_minimum_required(VERSION 3.24)
project(kache_hash_consumer LANGUAGES C CXX)

set(KACHE_HASH_SOURCE "" CACHE PATH "kache-hash's source tree, to take through FetchContent rather than find_package.")
if(KACHE_HASH_SOURCE)
    include(FetchContent)
    FetchContent_Declare(kache_hash SOURCE_DIR ${KACHE_HASH_SOURCE} FIND_PACKAGE_ARGS COMPONENTS capi counter)
    FetchContent_MakeAvailable(kache_hash)
else()
    find_package(kache_hash REQUIRED COMPONENTS capi counter)
endif()

# The table, whose ISA flags are its consumer's.
add_executable(table table.cpp)
target_compile_options(table PRIVATE -mavx2)
target_link_libraries(table PRIVATE kache_hash::kache_hash)

# The C interface, from C99.
add_executable(capi capi.c)
set_target_properties(capi PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED ON)
target_link_libraries(capi PRIVATE kache_hash::kache_capi)

# The counter, from C++17, with no AVX2.
add_executable(counter counter.cpp)
set_target_properties(counter PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)
target_link_libraries(counter PRIVATE kache_hash::kache_counter)
