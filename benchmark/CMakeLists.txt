
# Include paths.
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)


# Executables.
add_executable(measure measure.cpp)

# Benchmark executables, one per table.
add_executable(bench_kache_hash bench-kache-hash.cpp)


# Configs for targets.

target_include_directories(measure PRIVATE ${INCLUDE_DIR})
target_include_directories(bench_kache_hash PRIVATE ${INCLUDE_DIR})

target_compile_options(measure PRIVATE ${COMPILE_FLAGS})
target_compile_options(bench_kache_hash PRIVATE ${COMPILE_FLAGS})

target_link_libraries(measure PRIVATE ${LINK_FLAGS})
target_link_libraries(bench_kache_hash PRIVATE ${LINK_FLAGS})


# Target-links.

target_link_libraries(measure PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads)
target_link_libraries(bench_kache_hash PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads)


# Installation path.
install(TARGETS measure bench_kache_hash
        DESTINATION bin
)
