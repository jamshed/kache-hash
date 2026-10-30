
# Minimum required version of CMake.
cmake_minimum_required(VERSION 3.14)

# Project name, its version, and languages used.
set(PROJECT_NAME kache-hash)
project(${PROJECT_NAME}
            VERSION 0.5.0
            LANGUAGES CXX C
        )

# Whether this is the top-level project, and not a subproject of one that links the `kache_hash` target, e.g. through
# `FetchContent`: only top-level builds set the build type, export compile commands, and build the benchmarks.
string(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}" KACHE_HASH_TOP_LEVEL)

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

# Optimization flags: https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html. An x86-64 build on Apple silicon, to run
# under Rosetta, can't target the native CPU, and targets Haswell's AVX2 instead.
if(APPLE AND CMAKE_OSX_ARCHITECTURES STREQUAL "x86_64" AND CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "arm64")
    set(OPTIMIZE_FLAGS -funroll-loops -march=haswell)
else()
    set(OPTIMIZE_FLAGS -funroll-loops -march=native)
endif()

# Keep branches from crossing or ending on 32-byte boundaries: Skylake-derived cores (such as Cascade Lake) can't serve
# those from the decoded-uop cache under their JCC-erratum microcode. Other cores pay only for the padding. On Cascade
# Lake it gains 1-4% in builds and up to 9% in negative queries.
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|amd64")
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(JCC_PAD_FLAG -mbranches-within-32B-boundaries)
    else()
        set(JCC_PAD_FLAG -Wa,-mbranches-within-32B-boundaries)
    endif()

    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag(${JCC_PAD_FLAG} HAS_JCC_PAD)
    if(HAS_JCC_PAD)
        list(APPEND OPTIMIZE_FLAGS ${JCC_PAD_FLAG})
    endif()
endif()


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

if(AVX_512)
    add_compile_definitions(AVX_512)
endif()


# Required preprocessor definitions.
set(XXHASH_DEFS XXH_INLINE_ALL XXH_STATIC_LINKING_ONLY XXH_FORCE_MEMORY_ACCESS=0)   # `memcpy` reads: xxHash's default for GCC breaks strict aliasing.
add_compile_definitions(${XXHASH_DEFS})
add_compile_definitions(DEVELOP)
if(KACHE_HASH_TOP_LEVEL)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)   # Better `clangd` scan.
endif()

execute_process(
    COMMAND getconf LEVEL1_DCACHE_LINESIZE
    COMMAND tr -d '\n'
    OUTPUT_VARIABLE L1_CACHE_LINE_SIZE
    ERROR_QUIET
)

if(NOT L1_CACHE_LINE_SIZE MATCHES "^[1-9][0-9]*$")  # E.g. on macOS, whose `getconf` lacks the parameter.
    message("Cache-line size unknown; assuming 64 bytes.")
    set(L1_CACHE_LINE_SIZE 64)
endif()

add_compile_definitions(L1_CACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE})


# Required packages / libraries from the system.
find_package(Threads REQUIRED)


# The table as a header-only library, for other projects to link. It carries what its headers need: C++20, threads, the
# include path, and the xxHash and cache-line definitions. Optimization and ISA flags (AVX2 at least) are the linking
# project's own.
add_library(kache_hash INTERFACE)
add_library(kache_hash::kache_hash ALIAS kache_hash)
target_include_directories(kache_hash INTERFACE ${PROJECT_SOURCE_DIR}/include)
target_compile_definitions(kache_hash INTERFACE ${XXHASH_DEFS} L1_CACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE})
target_compile_features(kache_hash INTERFACE cxx_std_20)
target_link_libraries(kache_hash INTERFACE Threads::Threads)

# The k-mer counting library, in `counter/`: for the projects that link it, and in top-level builds if `COUNTER` is set.
if(COUNTER OR NOT KACHE_HASH_TOP_LEVEL)
    add_subdirectory(counter)
endif()


# The rest is for the benchmarks, in top-level builds only: in a subproject, it would set the parent's build type, write
# into the parent's source tree, and fetch the baselines.
if(NOT KACHE_HASH_TOP_LEVEL)
    return()
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(ZLIB REQUIRED)

# Required external projects.
include(ExternalProject)

set(EXT_LIB ${CMAKE_SOURCE_DIR}/external/lib/)
set(EXT_INCLUDE ${CMAKE_SOURCE_DIR}/external/include/)

file(MAKE_DIRECTORY ${EXT_LIB} ${EXT_INCLUDE})


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


# Build type.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build Type" FORCE)
endif()


add_subdirectory(benchmark)
