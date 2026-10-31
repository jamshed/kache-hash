
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


# Prepare the `Jellyfish` library. It provides a concurrent k-mer counting table. Its own build system configures and
# builds it, in the build tree: its headers take their fast paths, e.g. for hashing, only with the `config.h` that its
# `configure` writes there.
set(JELLYFISH_VERSION "2.3.1")
set(JELLYFISH_DIR ${CMAKE_SOURCE_DIR}/external/jellyfish-${JELLYFISH_VERSION})
set(JELLYFISH_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/jellyfish)
message("Build system will fetch and build the `Jellyfish` library")
ExternalProject_Add(prj_jellyfish
    DOWNLOAD_DIR        ${CMAKE_SOURCE_DIR}/external
    DOWNLOAD_COMMAND    curl -k -L https://github.com/gmarcais/Jellyfish/releases/download/v${JELLYFISH_VERSION}/jellyfish-${JELLYFISH_VERSION}.tar.gz -o jellyfish-${JELLYFISH_VERSION}.tar.gz &&
                        tar -xzf jellyfish-${JELLYFISH_VERSION}.tar.gz &&
                        rm jellyfish-${JELLYFISH_VERSION}.tar.gz

    SOURCE_DIR          ${JELLYFISH_DIR}
    BINARY_DIR          ${JELLYFISH_BUILD_DIR}
    CONFIGURE_COMMAND   ${JELLYFISH_DIR}/configure CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
                        --disable-shared --silent
    BUILD_COMMAND       $(MAKE) libjellyfish-2.0.la   # Its library only, not its tools.
    INSTALL_COMMAND     cp .libs/libjellyfish-2.0.a ${EXT_LIB}
)

add_library(jellyfish STATIC IMPORTED)
set_target_properties(jellyfish PROPERTIES IMPORTED_LOCATION ${EXT_LIB}/libjellyfish-2.0.a)
add_dependencies(jellyfish prj_jellyfish)


# Executables.
add_executable(measure measure.cpp)

# Benchmark executables, one per table.
add_executable(bench_kache_hash bench-kache-hash.cpp)
add_executable(bench_boost bench-boost.cpp)
add_executable(bench_libcuckoo bench-libcuckoo.cpp)
add_executable(bench_jellyfish bench-jellyfish.cpp)


# Dependencies.

add_dependencies(bench_boost prj_boost)
add_dependencies(bench_libcuckoo prj_libcuckoo)


# Configs for targets.

target_include_directories(measure PRIVATE ${INCLUDE_DIR})
target_include_directories(bench_kache_hash PRIVATE ${INCLUDE_DIR})
target_include_directories(bench_boost PRIVATE ${INCLUDE_DIR})
target_include_directories(bench_libcuckoo PRIVATE ${INCLUDE_DIR})
target_include_directories(bench_jellyfish PRIVATE ${INCLUDE_DIR})
target_include_directories(bench_boost SYSTEM PRIVATE ${EXT_INCLUDE})   # Third-party headers; their warnings are not ours.
target_include_directories(bench_libcuckoo SYSTEM PRIVATE ${EXT_INCLUDE})
target_include_directories(bench_jellyfish SYSTEM PRIVATE ${JELLYFISH_DIR}/include ${JELLYFISH_BUILD_DIR})

target_compile_definitions(bench_jellyfish PRIVATE HAVE_CONFIG_H)  # Jellyfish's headers then include its `config.h`.

target_compile_options(measure PRIVATE ${COMPILE_FLAGS})
target_compile_options(bench_kache_hash PRIVATE ${COMPILE_FLAGS})
target_compile_options(bench_boost PRIVATE ${COMPILE_FLAGS})
target_compile_options(bench_libcuckoo PRIVATE ${COMPILE_FLAGS})
target_compile_options(bench_jellyfish PRIVATE ${COMPILE_FLAGS})

target_link_libraries(measure PRIVATE ${LINK_FLAGS})
target_link_libraries(bench_kache_hash PRIVATE ${LINK_FLAGS})
target_link_libraries(bench_boost PRIVATE ${LINK_FLAGS})
target_link_libraries(bench_libcuckoo PRIVATE ${LINK_FLAGS})
target_link_libraries(bench_jellyfish PRIVATE ${LINK_FLAGS})


# Target-links.

target_link_libraries(measure PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads)
target_link_libraries(bench_kache_hash PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads)
target_link_libraries(bench_boost PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads)
target_link_libraries(bench_libcuckoo PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads)
target_link_libraries(bench_jellyfish PRIVATE jellyfish jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads)


# Installation path.
install(TARGETS measure bench_kache_hash bench_boost bench_libcuckoo bench_jellyfish
        DESTINATION bin
)


# kache-hash's benchmark at each minimizer size l in 15–31, as `bench_kache_hash_l<l>`, if `LSWEEP` is set: the drivers
# of the minimizer-size sweep.
if(LSWEEP)
    foreach(L RANGE 15 31)
        add_executable(bench_kache_hash_l${L} bench-kache-hash.cpp)

        target_compile_definitions(bench_kache_hash_l${L} PRIVATE MINIMIZER_SIZE=${L})
        target_include_directories(bench_kache_hash_l${L} PRIVATE ${INCLUDE_DIR})
        target_compile_options(bench_kache_hash_l${L} PRIVATE ${COMPILE_FLAGS})
        target_link_libraries(bench_kache_hash_l${L} PRIVATE ${LINK_FLAGS})
        target_link_libraries(bench_kache_hash_l${L} PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads)

        install(TARGETS bench_kache_hash_l${L}
                DESTINATION bin
        )
    endforeach()
endif()


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


# Rust's package manager, for the Rust projects of the baselines.
if(CBL OR SSHASH)
    find_program(CARGO cargo HINTS $ENV{CARGO_HOME}/bin $ENV{HOME}/.cargo/bin)
    if(NOT CARGO)
        message(FATAL_ERROR "Building the CBL and SSHash benchmarks requires `cargo`, from `rustup`.")
    endif()
endif()


# The CBL set's benchmark, if `CBL` is set. CBL is a Rust library; the crate in `cbl/` builds it, with a C interface,
# into a static library. Building it needs `rustup`, for the nightly toolchain pinned in `cbl/rust-toolchain.toml`, and
# libclang, which is looked up in `LIBCLANG_PATH` first. A libclang without its builtin headers, as anaconda's, needs
# some through `BINDGEN_EXTRA_CLANG_ARGS`, e.g. `-isystem $(gcc -print-file-name=include)`.
if(CBL)
    message("Build system will build the CBL library")
    ExternalProject_Add(prj_cbl
        SOURCE_DIR          ${CMAKE_CURRENT_SOURCE_DIR}/cbl
        BUILD_IN_SOURCE     TRUE    # For `rustup` to pick the pinned toolchain.
        CONFIGURE_COMMAND   ""
        # Make's jobserver is not passed down; Cargo warns of it in `MAKEFLAGS` and `MFLAGS`.
        BUILD_COMMAND       ${CMAKE_COMMAND} -E env --unset=MAKEFLAGS --unset=MFLAGS RUSTFLAGS=-Ctarget-cpu=native
                            ${CARGO} build --release --locked --target-dir ${CMAKE_CURRENT_BINARY_DIR}/cbl
        BUILD_ALWAYS        TRUE    # Cargo tracks the crate's changes itself.
        INSTALL_COMMAND     ${CMAKE_COMMAND} -E copy_if_different
                            ${CMAKE_CURRENT_BINARY_DIR}/cbl/release/libbench_cbl.a ${EXT_LIB}
    )

    add_library(cbl STATIC IMPORTED)
    set_target_properties(cbl PROPERTIES IMPORTED_LOCATION ${EXT_LIB}/libbench_cbl.a)
    add_dependencies(cbl prj_cbl)

    add_executable(bench_cbl bench-cbl.cpp)

    target_include_directories(bench_cbl PRIVATE ${INCLUDE_DIR})
    target_compile_options(bench_cbl PRIVATE ${COMPILE_FLAGS})
    target_link_libraries(bench_cbl PRIVATE ${LINK_FLAGS})
    target_link_libraries(bench_cbl PRIVATE cbl jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads)

    install(TARGETS bench_cbl
            DESTINATION bin
    )
endif()


# The SSHash dictionary's benchmark, if `SSHASH` is set. SSHash is static: `scripts/bench/sshash-index.sh` builds an
# index with the `ggcat` and `sshash` tools built here, and `bench_sshash` loads it.
if(SSHASH)
    # GGCAT builds the duplicate-free input of SSHash: the maximal unitigs of the k-mers.
    set(GGCAT_VERSION "2.2.0")

    message("Build system will fetch and build GGCAT")
    ExternalProject_Add(prj_ggcat
        DOWNLOAD_DIR        ${CMAKE_SOURCE_DIR}/external
        DOWNLOAD_COMMAND    git clone --recursive --depth 1 --branch v${GGCAT_VERSION}
                            https://github.com/algbio/ggcat.git ggcat-${GGCAT_VERSION}
        SOURCE_DIR          ${CMAKE_SOURCE_DIR}/external/ggcat-${GGCAT_VERSION}
        BUILD_IN_SOURCE     TRUE    # For `rustup` to pick its toolchain.
        CONFIGURE_COMMAND   ""
        BUILD_COMMAND       ${CMAKE_COMMAND} -E env --unset=MAKEFLAGS --unset=MFLAGS RUSTFLAGS=-Ctarget-cpu=native
                            ${CARGO} build --release --locked --package ggcat_cmdline
                            --target-dir ${CMAKE_CURRENT_BINARY_DIR}/ggcat-target
        INSTALL_COMMAND     cp ${CMAKE_CURRENT_BINARY_DIR}/ggcat-target/release/ggcat ${CMAKE_CURRENT_BINARY_DIR}
    )

    set(SSHASH_VERSION "6.0.0")
    set(SSHASH_DIR ${CMAKE_SOURCE_DIR}/external/sshash-${SSHASH_VERSION})

    message("Build system will fetch and build SSHash")
    ExternalProject_Add(prj_sshash
        DOWNLOAD_DIR        ${CMAKE_SOURCE_DIR}/external
        DOWNLOAD_COMMAND    git clone --recursive --depth 1 --branch v${SSHASH_VERSION}
                            https://github.com/jermp/sshash.git sshash-${SSHASH_VERSION}
        SOURCE_DIR          ${SSHASH_DIR}
        CMAKE_ARGS          -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER} -D CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                            -D CMAKE_BUILD_TYPE=Release
                            -D SSHASH_USE_ARCH_NATIVE=ON
        BUILD_COMMAND       ${CMAKE_COMMAND} --build . --target sshash --parallel   # Its library, and its tool.
        INSTALL_COMMAND     cp libsshash_static.a ${EXT_LIB} && cp sshash ${CMAKE_CURRENT_BINARY_DIR}
    )

    add_library(sshash_static STATIC IMPORTED)
    set_target_properties(sshash_static PROPERTIES IMPORTED_LOCATION ${EXT_LIB}/libsshash_static.a)
    add_dependencies(sshash_static prj_sshash)

    add_executable(bench_sshash bench-sshash.cpp)

    target_include_directories(bench_sshash PRIVATE ${INCLUDE_DIR})
    target_include_directories(bench_sshash SYSTEM PRIVATE  # As SSHash's own targets do.
                                ${SSHASH_DIR} ${SSHASH_DIR}/include
                                ${SSHASH_DIR}/external/pthash/include
                                ${SSHASH_DIR}/external/pthash/external/bits/include
                                ${SSHASH_DIR}/external/pthash/external/bits/external/essentials/include
                                ${SSHASH_DIR}/external/pthash/external/fastmod
                                ${SSHASH_DIR}/external/pthash/external/xxHash
                                ${SSHASH_DIR}/external/pthash/external/mm_file/include)
    target_compile_options(bench_sshash PRIVATE ${COMPILE_FLAGS})
    target_link_libraries(bench_sshash PRIVATE ${LINK_FLAGS})
    target_link_libraries(bench_sshash PRIVATE sshash_static jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads)

    install(TARGETS bench_sshash
            DESTINATION bin
    )
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ggcat ${CMAKE_CURRENT_BINARY_DIR}/sshash
            DESTINATION bin
    )
endif()
