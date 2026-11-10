
# Builds kache-hash's C interface library for the crate: kache-hash as a subproject, so that only its libraries are
# configured, with none of its benchmarks; and only the C interface library is built, then installed for `build.rs`.

cmake_minimum_required(VERSION 3.14)
project(kache_hash_sys LANGUAGES C CXX)

set(KACHE_HASH_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../.." CACHE PATH "Root of the kache-hash source tree")

# The C interface library alone, and static: the crate links it into the Rust binary.
set(BUILD_SHARED_LIBS OFF)
set(KACHE_HASH_BUILD_COUNTER OFF)

# The library's targets are built only on demand: `all` builds the C interface alone.
add_subdirectory(${KACHE_HASH_ROOT} kache-hash EXCLUDE_FROM_ALL)
add_custom_target(kache_capi_all ALL)
add_dependencies(kache_capi_all kache_capi)

install(TARGETS kache_capi ARCHIVE DESTINATION lib)
