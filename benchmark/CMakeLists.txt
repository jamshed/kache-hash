
# Include paths.
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)


# Baselines' external projects.

# Prepare the `boost` library. It provides scalable concurrent hash tables. The version is pinned, for benchmarks to
# be reproducible across machines.
set(BOOST_MIN_VERSION "1.89.0")
message("Build system will fetch and install required Boost libraries")
ExternalProject_Add(prj_boost
    DOWNLOAD_DIR        ${CMAKE_SOURCE_DIR}/external
    DOWNLOAD_COMMAND    curl -k -L https://github.com/boostorg/boost/releases/download/boost-${BOOST_MIN_VERSION}/boost-${BOOST_MIN_VERSION}-cmake.tar.gz -o boost-${BOOST_MIN_VERSION}-cmake.tar.gz &&
                        tar -xzf boost-${BOOST_MIN_VERSION}-cmake.tar.gz &&
                        rm boost-${BOOST_MIN_VERSION}-cmake.tar.gz

    SOURCE_DIR          ${CMAKE_SOURCE_DIR}/external/boost-${BOOST_MIN_VERSION}
    BUILD_IN_SOURCE     TRUE
    INSTALL_DIR         ${CMAKE_SOURCE_DIR}/external
    CONFIGURE_COMMAND   mkdir -p build && cd build &&
                        cmake -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER} -D CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                        -D BOOST_INCLUDE_LIBRARIES=unordered
                        -D CMAKE_BUILD_TYPE=Release
                        -D CMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/boost-${BOOST_MIN_VERSION}/build/ ..
    BUILD_COMMAND       cd build && cmake --build . --target install --config Release
    INSTALL_COMMAND     cp -rf build/include/boost ${EXT_INCLUDE}/
)


# Prepare the `libcuckoo` library.
set(LIBCUCKOO_MIN_VERSION "0.3.1")
message("Build system will fetch and install the `libcuckoo` library")
ExternalProject_Add(prj_libcuckoo
    DOWNLOAD_DIR        ${CMAKE_SOURCE_DIR}/external
    DOWNLOAD_COMMAND    curl -k -L https://github.com/efficient/libcuckoo/archive/refs/tags/v${LIBCUCKOO_MIN_VERSION}.tar.gz -o libcuckoo-${LIBCUCKOO_MIN_VERSION}.tar.gz &&
                        tar -xzf libcuckoo-${LIBCUCKOO_MIN_VERSION}.tar.gz &&
                        rm libcuckoo-${LIBCUCKOO_MIN_VERSION}.tar.gz

    SOURCE_DIR          ${CMAKE_SOURCE_DIR}/external/libcuckoo-${LIBCUCKOO_MIN_VERSION}
    BUILD_IN_SOURCE     TRUE
    INSTALL_DIR         ${CMAKE_SOURCE_DIR}/external
    CONFIGURE_COMMAND   ""
    BUILD_COMMAND       ""
    INSTALL_COMMAND     cp -rf libcuckoo ${EXT_INCLUDE}/
)


# Executables.
add_executable(measure measure.cpp)

# Benchmark executables, one per table.
add_executable(bench_kache_hash bench-kache-hash.cpp)
add_executable(bench_boost bench-boost.cpp)
add_executable(bench_libcuckoo bench-libcuckoo.cpp)


# Dependencies.

add_dependencies(bench_boost prj_boost)
add_dependencies(bench_libcuckoo prj_libcuckoo)


# Configs for targets.

target_include_directories(measure PRIVATE ${INCLUDE_DIR})
target_include_directories(bench_kache_hash PRIVATE ${INCLUDE_DIR})
target_include_directories(bench_boost PRIVATE ${INCLUDE_DIR})
target_include_directories(bench_libcuckoo PRIVATE ${INCLUDE_DIR})
target_include_directories(bench_boost SYSTEM PRIVATE ${EXT_INCLUDE})   # Third-party headers; their warnings are not ours.
target_include_directories(bench_libcuckoo SYSTEM PRIVATE ${EXT_INCLUDE})

target_compile_options(measure PRIVATE ${COMPILE_FLAGS})
target_compile_options(bench_kache_hash PRIVATE ${COMPILE_FLAGS})
target_compile_options(bench_boost PRIVATE ${COMPILE_FLAGS})
target_compile_options(bench_libcuckoo PRIVATE ${COMPILE_FLAGS})

target_link_libraries(measure PRIVATE ${LINK_FLAGS})
target_link_libraries(bench_kache_hash PRIVATE ${LINK_FLAGS})
target_link_libraries(bench_boost PRIVATE ${LINK_FLAGS})
target_link_libraries(bench_libcuckoo PRIVATE ${LINK_FLAGS})


# Target-links.

target_link_libraries(measure PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads)
target_link_libraries(bench_kache_hash PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads)
target_link_libraries(bench_boost PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads)
target_link_libraries(bench_libcuckoo PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads)


# Installation path.
install(TARGETS measure bench_kache_hash bench_boost bench_libcuckoo
        DESTINATION bin
)


# The iceberg table's benchmark, on Linux only. It builds the vendored iceberg sources, in `icebergHT/`, as iceberg's
# own Makefile does by default—vectorized, with block locks and resizes—and links them with the driver under LTO.
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    file(GLOB ICEBERG_SRC ${CMAKE_CURRENT_SOURCE_DIR}/icebergHT/*.c)
    set_source_files_properties(${ICEBERG_SRC} PROPERTIES COMPILE_OPTIONS -w)  # Third-party sources; their warnings are not ours.

    add_executable(bench_iceberg bench-iceberg.cpp ${ICEBERG_SRC})

    target_compile_definitions(bench_iceberg PRIVATE USE_AVX ENABLE_BLOCK_LOCKING ENABLE_RESIZE)
    target_include_directories(bench_iceberg PRIVATE ${INCLUDE_DIR})
    target_compile_options(bench_iceberg PRIVATE ${COMPILE_FLAGS})
    target_link_libraries(bench_iceberg PRIVATE ${LINK_FLAGS})
    target_link_libraries(bench_iceberg PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads)

    include(CheckIPOSupported)
    check_ipo_supported(RESULT IPO_SUPPORTED OUTPUT IPO_ERROR)
    if(IPO_SUPPORTED)
        set_property(TARGET bench_iceberg PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    else()
        message(WARNING "Building `bench_iceberg` without LTO: ${IPO_ERROR}")
    endif()

    install(TARGETS bench_iceberg
            DESTINATION bin
    )
endif()
