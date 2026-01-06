
# Minimum required version of CMake.
cmake_minimum_required(VERSION 3.14)

# Project name, its version, and languages used.
set(PROJECT_NAME kache-hash)
project(${PROJECT_NAME}
            VERSION 0.1.0
            LANGUAGES CXX C
        )

# Language standards.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Platform-specific compiler standards.
set(CMAKE_CXX_EXTENSIONS OFF)

# Compilers for sanitizer-support.
if(ASAN OR UBSAN)
    set(CLANG TRUE)
endif()

if(CLANG)
    set(CMAKE_C_COMPILER clang)
    set(CMAKE_CXX_COMPILER clang++)
    # set(CMAKE_LINKER_TYPE LLD)
endif()

# Minimum compiler versions.
set(GCC_VERSION_MIN "9.1")
set(CLANG_VERSION_MIN "9.0")

if(ASAN OR UBSAN)
    set(CLANG TRUE)
endif()

if(CLANG)
    set(CMAKE_C_COMPILER clang)
    set(CMAKE_CXX_COMPILER clang++)
    # set(CMAKE_LINKER_TYPE LLD)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS GCC_VERSION_MIN)
        message(FATAL_ERROR "${PROJECT_NAME} requires GCC version to be at least ${GCC_VERSION_MIN}."
                            "Available version is ${CMAKE_CXX_COMPILER_VERSION}.")
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS CLANG_VERSION_MIN)
        message(FATAL_ERROR "${PROJECT_NAME} requires Clang version to be at least ${CLANG_VERSION_MIN}."
                            "Available version is ${CMAKE_CXX_COMPILER_VERSION}.")
    endif()
endif()


# Warning flags: https://gcc.gnu.org/onlinedocs/gcc/Warning-Options.html.
set(WARNING_FLAGS   -Wall -Wextra -Wpedantic)
set(SUPPRESS_WARNING_FLAGS)

# Warning flags to be suppressed with clang.
set(CLANG_SUPPRESS_WARNING_FLAGS    -Wno-unknown-attributes)

if(CLANG)
    list(APPEND SUPPRESS_WARNING_FLAGS ${CLANG_SUPPRESS_WARNING_FLAGS})
endif()

# Debug flags: https://gcc.gnu.org/onlinedocs/gcc/Debugging-Options.html.
set(DEBUG_FLAGS -ggdb)

# Address sanitizer flags: https://github.com/google/sanitizers/wiki/AddressSanitizer.
set(ASAN_FLAGS -fsanitize=address -fno-omit-frame-pointer)
# set(ASAN_LINK_FLAGS -static-libasan)

# Undefined behavior sanitizer flags: https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html.
set(UBSAN_FLAGS -fsanitize=undefined
                -fsanitize=float-divide-by-zero #-fsanitize=unsigned-integer-overflow
                -fsanitize=implicit-conversion -fsanitize=local-bounds -fsanitize=nullability)
# set(UBSAN_LINK_FLAGS -fuse-ld=lld)

# Optimization flags: https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html.
set(OPTIMIZE_FLAGS -funroll-loops -march=native)


# Compile- and link-flags.
set(COMPILE_FLAGS "")
set(LINK_FLAGS "")

list(APPEND COMPILE_FLAGS ${WARNING_FLAGS} ${SUPPRESS_WARNING_FLAGS} ${OPTIMIZE_FLAGS})

if(DEBUG_INFO)
    list(APPEND COMPILE_FLAGS ${DEBUG_FLAGS})
endif()

if(ASAN)
    list(APPEND COMPILE_FLAGS ${ASAN_FLAGS} ${DEBUG_FLAGS})
    list(APPEND LINK_FLAGS ${ASAN_LINK_FLAGS} ${ASAN_FLAGS})
endif()

if(UBSAN)
    list(APPEND COMPILE_FLAGS ${UBSAN_FLAGS} ${DEBUG_FLAGS})
    list(APPEND LINK_FLAGS ${UBSAN_LINK_FLAGS} ${UBSAN_FLAGS})
endif()

if(QUIET)
    add_compile_definitions(QUIET)
endif()

if(AVX_512)
    add_compile_definitions(AVX_512)
endif()


# Required preprocessor definitions.
set(XXHASH_DEFS XXH_INLINE_ALL XXH_STATIC_LINKING_ONLY)
add_compile_definitions(${XXHASH_DEFS})
add_compile_definitions(DEVELOP)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)   # Better `clangd` scan.

execute_process(
    COMMAND getconf LEVEL1_DCACHE_LINESIZE
    COMMAND tr -d '\n'
    OUTPUT_VARIABLE L1_CACHE_LINE_SIZE
)

add_compile_definitions(L1_CACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE})


# Required packages / libraries from the system.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)
# find_package(PFM REQUIRED)

# Required external projects.
include(ExternalProject)

set(EXT_LIB ${CMAKE_SOURCE_DIR}/external/lib/)
set(EXT_INCLUDE ${CMAKE_SOURCE_DIR}/external/include/)

file(MAKE_DIRECTORY ${EXT_LIB} ${EXT_INCLUDE})

#[[
set(BENCHMARK_MIN_VERSION "v1.9.4")
message("Build system will fetch and install Google Benchmark.")
ExternalProject_Add(prj_benchmark
    DOWNLOAD_DIR        ${CMAKE_SOURCE_DIR}/external
    DOWNLOAD_COMMAND    curl -k -L https://github.com/google/benchmark/archive/refs/tags/${BENCHMARK_MIN_VERSION}.tar.gz -o benchmark-${BENCHMARK_MIN_VERSION}.tar.gz &&
                        tar -xzf benchmark-${BENCHMARK_MIN_VERSION}.tar.gz -C benchmark-${BENCHMARK_MIN_VERSION} --strip-components 1 &&
                        rm benchmark-${BENCHMARK_MIN_VERSION}.tar.gz

    SOURCE_DIR          ${CMAKE_SOURCE_DIR}/external/benchmark-${BENCHMARK_MIN_VERSION}
    BUILD_IN_SOURCE     TRUE
    INSTALL_DIR         ${CMAKE_SOURCE_DIR}/external
    CONFIGURE_COMMAND   cmake -E make_directory "build" &&
                        cmake -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER} -D CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                        -D BENCHMARK_DOWNLOAD_DEPENDENCIES=on -D CMAKE_BUILD_TYPE=Release -D BENCHMARK_ENABLE_LIBPFM=ON
                        -S . -B "build"
    BUILD_COMMAND       cmake --build "build" --config Release
    INSTALL_COMMAND     cp build/src/libbenchmark.a ${EXT_LIB} &&
                        cp -rf include/benchmark ${EXT_INCLUDE}
)

set(BENCHMARK_LIB ${EXT_LIB}/libbenchmark.a)

add_library(benchmark STATIC IMPORTED)
set_target_properties(benchmark PROPERTIES IMPORTED_LOCATION ${BENCHMARK_LIB})
add_dependencies(benchmark prj_benchmark)
#]]


set(MALLOC_LIB "")
set(JEMALLOC_MIN_VERSION "5.3.0")
find_package(Jemalloc ${JEMALLOC_MIN_VERSION})
if(JEMALLOC_FOUND)
    message("Found jemalloc (v${JEMALLOC_VERSION}) in the system.")
    set(MALLOC_LIB ${JEMALLOC_LIBRARIES})
else()
    message("Build system will fetch and install jemalloc.")
    ExternalProject_Add(prj_jemalloc
        DOWNLOAD_DIR        ${CMAKE_SOURCE_DIR}/external
        DOWNLOAD_COMMAND    curl -k -L https://github.com/jemalloc/jemalloc/archive/${JEMALLOC_MIN_VERSION}.tar.gz -o jemalloc-${JEMALLOC_MIN_VERSION}.tar.gz &&
                            tar -xzf jemalloc-${JEMALLOC_MIN_VERSION}.tar.gz &&
                            rm jemalloc-${JEMALLOC_MIN_VERSION}.tar.gz

        SOURCE_DIR          ${CMAKE_SOURCE_DIR}/external/jemalloc-${JEMALLOC_MIN_VERSION}
        BUILD_IN_SOURCE     TRUE
        INSTALL_DIR         ${CMAKE_SOURCE_DIR}/external
        CONFIGURE_COMMAND   sh -c "CC=${CMAKE_C_COMPILER} ./autogen.sh --disable-debug --enable-static --prefix=<INSTALL_DIR> --silent"
        INSTALL_COMMAND     cp lib/libjemalloc.a ${EXT_LIB}
    )

    set(MALLOC_LIB ${EXT_LIB}/libjemalloc.a)
endif()

add_library(jemalloc STATIC IMPORTED)
set_target_properties(jemalloc PROPERTIES IMPORTED_LOCATION ${MALLOC_LIB})
if(NOT JEMALLOC_FOUND)
    add_dependencies(jemalloc prj_jemalloc)
endif()


# Prepare the `boost` library. It provides a scalable concurrent hashmap.
set(BOOST_MIN_VERSION "1.89.0")
find_package(Boost ${BOOST_MIN_VERSION} CONFIG)
if(BOOST_FOUND)
    message("Found boost (v${BOOST_VERSION}) in the system at ${BOOST_LIBRARIES}")
else()
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
endif()


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
    INSTALL_COMMAND     cp -rf libcuckoo/ ${EXT_INCLUDE}/
)


# Prepare the `parlay` library—a parallel algorithms programming toolkit.
message("Build system will fetch and install parlaylib")
ExternalProject_Add(prj_parlaylib
    DOWNLOAD_DIR        ${CMAKE_SOURCE_DIR}/external
    DOWNLOAD_COMMAND    git clone https://github.com/cmuparlay/parlaylib.git
    SOURCE_DIR          ${CMAKE_SOURCE_DIR}/external/parlaylib
    BUILD_IN_SOURCE     TRUE
    INSTALL_DIR         ${CMAKE_SOURCE_DIR}/external
    CONFIGURE_COMMAND   ""
    BUILD_COMMAND       ""
    INSTALL_COMMAND     cp -rf include/parlay ${EXT_INCLUDE}
)


#[[
set(NTHASH_LIB "")
set(NTHASH_MIN_VERSION "v2.3.0")
message("Build system will fetch and install ntHash for benchmarking.")
ExternalProject_Add(prj_nthash
    DOWNLOAD_DIR        ${CMAKE_SOURCE_DIR}/external
    DOWNLOAD_COMMAND    git clone git@github.com:bcgsc/ntHash.git

    SOURCE_DIR          ${CMAKE_SOURCE_DIR}/external/ntHash
    BUILD_IN_SOURCE     TRUE
    INSTALL_DIR         ${CMAKE_SOURCE_DIR}/external
    CONFIGURE_COMMAND   meson setup --buildtype=release --prefix=${CMAKE_SOURCE_DIR}/external/ build
    BUILD_COMMAND       meson install -C build
    INSTALL_COMMAND     ""
)

set(NTHASH_LIB ${EXT_LIB}/libnthash.a)
add_library(nthash STATIC IMPORTED)
set_target_properties(nthash PROPERTIES IMPORTED_LOCATION ${NTHASH_LIB})
add_dependencies(nthash prj_nthash)
#]]


# Build type.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build Type" FORCE)
endif()


add_subdirectory(src)
