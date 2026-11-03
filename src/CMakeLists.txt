
# Include paths.
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)


# Targets.


# Libraries.

file(GLOB ICEBERG_HT_SRC
        "icebergHT/*"
)

add_library(iceberg_ht STATIC ${ICEBERG_HT_SRC})
set(ICEBERG_INCLUDE ${INCLUDE_DIR}/icebergHT)


# Executables.
add_executable(example examples/example.cpp)

# Benchmark experiment executables.
add_executable(throughput_sequential experiments/throughput-sequential.cpp)
add_executable(throughput_concurrent experiments/throughput-concurrent.cpp)
add_executable(cache_miss experiments/cache-miss.cpp)
add_executable(cache_miss_query experiments/cache-miss-query.cpp)
add_executable(query experiments/query.cpp)
add_executable(key_distribution experiments/key-distribution.cpp)
add_executable(minimizer_effect experiments/minimizer-effect.cpp)
add_executable(lf_v_key_dist experiments/lf-vs-key-dist.cpp)

# Test executables.
add_executable(concurrent_insert_test tests/concurrent-insert.cpp)


# Dependencies


# Configs for targets.

target_include_directories(example PRIVATE ${INCLUDE_DIR})
target_include_directories(iceberg_ht PRIVATE ${ICEBERG_INCLUDE})

target_include_directories(throughput_sequential PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})
target_include_directories(throughput_concurrent PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})
target_include_directories(cache_miss PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})
target_include_directories(cache_miss_query PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})
target_include_directories(query PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})
target_include_directories(key_distribution PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})
target_include_directories(minimizer_effect PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})
target_include_directories(lf_v_key_dist PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})
target_include_directories(concurrent_insert_test PRIVATE ${INCLUDE_DIR})

target_compile_options(example PRIVATE ${COMPILE_FLAGS})
target_compile_options(iceberg_ht PRIVATE ${COMPILE_FLAGS})

target_compile_options(throughput_sequential PRIVATE ${COMPILE_FLAGS})
target_compile_options(throughput_concurrent PRIVATE ${COMPILE_FLAGS})
target_compile_options(cache_miss PRIVATE ${COMPILE_FLAGS})
target_compile_options(cache_miss_query PRIVATE ${COMPILE_FLAGS})
target_compile_options(query PRIVATE ${COMPILE_FLAGS})
target_compile_options(key_distribution PRIVATE ${COMPILE_FLAGS})
target_compile_options(minimizer_effect PRIVATE ${COMPILE_FLAGS})
target_compile_options(lf_v_key_dist PRIVATE ${COMPILE_FLAGS})
target_compile_options(concurrent_insert_test PRIVATE ${COMPILE_FLAGS})

target_link_libraries(example PRIVATE ${LINK_FLAGS})
target_link_libraries(iceberg_ht PRIVATE ${LINK_FLAGS})

target_compile_options(throughput_sequential PRIVATE ${LINK_FLAGS})
target_compile_options(throughput_concurrent PRIVATE ${LINK_FLAGS})
target_compile_options(cache_miss PRIVATE ${LINK_FLAGS})
target_compile_options(cache_miss_query PRIVATE ${LINK_FLAGS})
target_compile_options(query PRIVATE ${LINK_FLAGS})
target_compile_options(key_distribution PRIVATE ${LINK_FLAGS})
target_compile_options(minimizer_effect PRIVATE ${LINK_FLAGS})
target_compile_options(lf_v_key_dist PRIVATE ${LINK_FLAGS})


# Target-links.

target_link_libraries(iceberg_ht PRIVATE jemalloc)
target_link_libraries(iceberg_ht PRIVATE ${CMAKE_DL_LIBS})    # For dynamic shared objects, like `jemalloc`.

target_link_libraries(example PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads)

target_link_libraries(throughput_sequential PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads iceberg_ht)
target_link_libraries(throughput_concurrent PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads iceberg_ht)
target_link_libraries(cache_miss PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads iceberg_ht)
target_link_libraries(cache_miss_query PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads iceberg_ht)
target_link_libraries(query PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads iceberg_ht)
target_link_libraries(key_distribution PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads iceberg_ht)
target_link_libraries(minimizer_effect PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads iceberg_ht)
target_link_libraries(lf_v_key_dist PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads iceberg_ht)
target_link_libraries(concurrent_insert_test PRIVATE ${LINK_FLAGS} ZLIB::ZLIB Threads::Threads)


# Installation path.
install(TARGETS example throughput_sequential throughput_concurrent cache_miss cache_miss_query query key_distribution minimizer_effect lf_v_key_dist
        DESTINATION bin
)
