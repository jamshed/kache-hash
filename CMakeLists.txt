
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

if(INSTRUMENT)
    add_compile_definitions(INSTRUMENT)
endif()

if(NO_TABLE_LOCK)
    add_compile_definitions(NO_TABLE_LOCK)
endif()


if(MIN_SIZE_EXP)
    add_compile_definitions(MIN_SIZE_EXP=${MIN_SIZE_EXP})
endif()

# Required preprocessor definitions.
set(XXHASH_DEFS XXH_INLINE_ALL XXH_STATIC_LINKING_ONLY XXH_FORCE_MEMORY_ACCESS=0)   # `memcpy` reads: xxHash's default for GCC breaks strict aliasing.
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


# Build type.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build Type" FORCE)
endif()


add_subdirectory(benchmark)
add_subdirectory(src)
