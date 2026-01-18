
# Include paths.
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)


# Targets.

# Libraries.

set(SRC
        Parser.cpp
        # examples/Kmer_Set.cpp
        # examples/Minimizer_Scan.cpp
)

file(GLOB ICEBERG_HT_SRC
        "icebergHT/*"
)

add_library(core STATIC ${SRC})

add_library(iceberg_ht STATIC ${ICEBERG_HT_SRC})
set(ICEBERG_INCLUDE ${INCLUDE_DIR}/icebergHT)


# Executables.
# add_executable(bench bench.cpp)
add_executable(bench benchmark.cpp)
add_executable(main main.cpp)

# Benchmark experiment executables.
add_executable(throughput_sequential experiments/throughput-sequential.cpp)
add_executable(throughput_concurrent experiments/throughput-concurrent.cpp)
add_executable(cache_miss experiments/cache-miss.cpp)
add_executable(cache_miss_query experiments/cache-miss-query.cpp)
add_executable(query experiments/query.cpp)
add_executable(key_distribution experiments/key-distribution.cpp)
add_executable(minimizer_effect experiments/minimizer-effect.cpp)


# Dependencies
add_dependencies(core prj_parlaylib)

# add_dependencies(core prj_benchmark)

if(NOT BOOST_FOUND)
    add_dependencies(core prj_boost)
endif()


# Configs for targets.

target_include_directories(core PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})
target_include_directories(iceberg_ht PRIVATE ${ICEBERG_INCLUDE})
target_include_directories(bench PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})
target_include_directories(main PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})

target_include_directories(throughput_sequential PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})
target_include_directories(throughput_concurrent PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})
target_include_directories(cache_miss PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})
target_include_directories(cache_miss_query PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})
target_include_directories(query PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})
target_include_directories(key_distribution PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})
target_include_directories(minimizer_effect PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})

target_compile_options(core PRIVATE ${COMPILE_FLAGS})
target_compile_options(iceberg_ht PRIVATE ${COMPILE_FLAGS})
target_compile_options(bench PRIVATE ${COMPILE_FLAGS})
target_compile_options(main PRIVATE ${COMPILE_FLAGS})

target_compile_options(throughput_sequential PRIVATE ${COMPILE_FLAGS})
target_compile_options(throughput_concurrent PRIVATE ${COMPILE_FLAGS})
target_compile_options(cache_miss PRIVATE ${COMPILE_FLAGS})
target_compile_options(cache_miss_query PRIVATE ${COMPILE_FLAGS})
target_compile_options(query PRIVATE ${COMPILE_FLAGS})
target_compile_options(key_distribution PRIVATE ${COMPILE_FLAGS})
target_compile_options(minimizer_effect PRIVATE ${COMPILE_FLAGS})

target_link_libraries(core PRIVATE ${LINK_FLAGS})
target_link_libraries(iceberg_ht PRIVATE ${LINK_FLAGS})
target_link_libraries(bench PRIVATE ${LINK_FLAGS})
target_link_libraries(main PRIVATE ${LINK_FLAGS})

target_compile_options(throughput_sequential PRIVATE ${LINK_FLAGS})
target_compile_options(throughput_concurrent PRIVATE ${LINK_FLAGS})
target_compile_options(cache_miss PRIVATE ${LINK_FLAGS})
target_compile_options(cache_miss_query PRIVATE ${LINK_FLAGS})
target_compile_options(query PRIVATE ${LINK_FLAGS})
target_compile_options(key_distribution PRIVATE ${LINK_FLAGS})
target_compile_options(minimizer_effect PRIVATE ${LINK_FLAGS})


# Target-links.
target_link_libraries(core PRIVATE jemalloc)
target_link_libraries(core PRIVATE ${CMAKE_DL_LIBS})    # For dynamic shared objects, like `jemalloc`.
target_link_libraries(core PRIVATE ZLIB::ZLIB)
target_link_libraries(core PRIVATE Threads::Threads)
# target_link_libraries(core PRIVATE benchmark)
# target_link_libraries(core PRIVATE pfm)

target_link_libraries(iceberg_ht PRIVATE jemalloc)
target_link_libraries(iceberg_ht PRIVATE ${CMAKE_DL_LIBS})    # For dynamic shared objects, like `jemalloc`.

target_link_libraries(bench PRIVATE core)
target_link_libraries(bench PRIVATE iceberg_ht)

target_link_libraries(throughput_sequential PRIVATE core iceberg_ht)
target_link_libraries(throughput_concurrent PRIVATE core iceberg_ht)
target_link_libraries(cache_miss PRIVATE core iceberg_ht)
target_link_libraries(cache_miss_query PRIVATE core iceberg_ht)
target_link_libraries(query PRIVATE core iceberg_ht)
target_link_libraries(key_distribution PRIVATE core iceberg_ht)
target_link_libraries(minimizer_effect PRIVATE core iceberg_ht)
# target_link_libraries(bench PRIVATE benchmark)
# target_link_libraries(bench PRIVATE nthash)

target_link_libraries(main PRIVATE core)
# target_link_libraries(main PRIVATE benchmark)
# target_link_libraries(main PRIVATE nthash)


# Installation path.
install(TARGETS main bench throughput_sequential throughput_concurrent cache_miss cache_miss_query query key_distribution minimizer_effect
        DESTINATION bin
        )
