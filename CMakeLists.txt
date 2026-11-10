
# Minimum required version of CMake.
cmake_minimum_required(VERSION 3.14)

# A parent project's normal variables keep their values over the cache entries that the options and settings below
# create, as they do over `option`s already.
if(POLICY CMP0126)
    cmake_policy(SET CMP0126 NEW)
endif()

# Project name, its version, and languages used. The version is the headers' `kache_hash/version.hpp`'s.
set(PROJECT_NAME kache-hash)
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/include/kache_hash/version.hpp" KACHE_HASH_VERSION_DEFS
     REGEX "^#define KACHE_HASH_VERSION_(MAJOR|MINOR|PATCH) [0-9]+$")
foreach(part MAJOR MINOR PATCH)
    if(NOT KACHE_HASH_VERSION_DEFS MATCHES "KACHE_HASH_VERSION_${part} ([0-9]+)")
        message(FATAL_ERROR "kache_hash/version.hpp lacks KACHE_HASH_VERSION_${part}.")
    endif()
    set(KACHE_HASH_VERSION_${part} ${CMAKE_MATCH_1})
endforeach()

project(${PROJECT_NAME}
            VERSION ${KACHE_HASH_VERSION_MAJOR}.${KACHE_HASH_VERSION_MINOR}.${KACHE_HASH_VERSION_PATCH}
            LANGUAGES CXX C
        )

# Whether this is the top-level project, and not a subproject of one that links its targets, e.g. through
# `FetchContent`: only top-level builds set the build type, export compile commands, install by default, and build the
# tests and the benchmarks.
string(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}" KACHE_HASH_TOP_LEVEL)


# Options.
option(KACHE_HASH_BUILD_CAPI "Build the C interface library, kache_capi." ON)
option(KACHE_HASH_BUILD_COUNTER "Build the k-mer counting library, kache_counter." ON)
option(KACHE_HASH_BUILD_BENCHMARKS "Build the benchmarks and the tools, fetching the baselines; top-level builds only." OFF)
option(KACHE_HASH_INSTALL "Generate the install rules." ${KACHE_HASH_TOP_LEVEL})
option(KACHE_HASH_BUILD_TESTS "Build the tests, for CTest; top-level builds only." ${KACHE_HASH_TOP_LEVEL})

# The ISA of the compiled targets, i.e. the libraries, the benchmarks and the tools: `native`, the building CPU's, by
# default; `x86-64-v3` (AVX2) or `x86-64-v4` (AVX-512), for binaries that run on other CPUs, with `x86-64-v3` the default
# under conda-build and cibuildwheel; or `none`, for the toolchain's flags alone. The header-only library takes its
# consumers' flags.
if(DEFINED ENV{CONDA_BUILD} OR DEFINED ENV{CIBUILDWHEEL})
    set(KACHE_HASH_ARCH_DEFAULT x86-64-v3)
else()
    set(KACHE_HASH_ARCH_DEFAULT native)
endif()
set(KACHE_HASH_ARCH ${KACHE_HASH_ARCH_DEFAULT} CACHE STRING "ISA of the compiled targets: native, x86-64-v3, x86-64-v4 or none.")
set_property(CACHE KACHE_HASH_ARCH PROPERTY STRINGS native x86-64-v3 x86-64-v4 none)

# Development options, for the compiled targets: `DEBUG_INFO`, and the address and undefined-behavior sanitizers, `ASAN`
# and `UBSAN`, the latter with Clang.

# The options of earlier releases, deprecated.
if(DEFINED KACHE_HASH_CAPI)
    message(DEPRECATION "KACHE_HASH_CAPI is deprecated: the C interface library builds unless KACHE_HASH_BUILD_CAPI is OFF.")
    set(KACHE_HASH_BUILD_CAPI ${KACHE_HASH_CAPI})
endif()
if(KACHE_HASH_TOP_LEVEL AND DEFINED COUNTER)
    message(DEPRECATION "COUNTER is deprecated: the counting library builds unless KACHE_HASH_BUILD_COUNTER is OFF.")
    set(KACHE_HASH_BUILD_COUNTER ${COUNTER})
endif()


# Release by default, in top-level builds.
if(KACHE_HASH_TOP_LEVEL AND NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type." FORCE)
endif()

# Language standards.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Platform-specific compiler standards.
set(CMAKE_CXX_EXTENSIONS OFF)

# Minimum compiler versions: GCC 11 and Clang 12, for C++20's `<barrier>` and the `-march` names of the x86-64 levels.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11)
    message(FATAL_ERROR "kache-hash needs GCC 11 or later; this is GCC ${CMAKE_CXX_COMPILER_VERSION}.")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 12)
    message(FATAL_ERROR "kache-hash needs Clang 12 or later; this is Clang ${CMAKE_CXX_COMPILER_VERSION}.")
endif()


# Whether the build targets x86-64: an x86-64 build on Apple silicon, to run under Rosetta, sets it through
# `CMAKE_OSX_ARCHITECTURES`.
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64)$" OR (APPLE AND CMAKE_OSX_ARCHITECTURES STREQUAL "x86_64"))
    set(KACHE_HASH_X86 TRUE)
else()
    set(KACHE_HASH_X86 FALSE)
endif()

# The ISA flags of `KACHE_HASH_ARCH`. Under Rosetta, which can't run the native CPU's code, `native` targets Haswell's
# AVX2.
if(KACHE_HASH_ARCH STREQUAL "native")
    if(APPLE AND CMAKE_OSX_ARCHITECTURES STREQUAL "x86_64" AND CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "arm64")
        set(ARCH_FLAGS -march=haswell)
    else()
        set(ARCH_FLAGS -march=native)
    endif()
elseif(KACHE_HASH_ARCH MATCHES "^x86-64-v[34]$")
    set(ARCH_FLAGS -march=${KACHE_HASH_ARCH})
elseif(KACHE_HASH_ARCH STREQUAL "none")
    set(ARCH_FLAGS "")
else()
    message(FATAL_ERROR "KACHE_HASH_ARCH is ${KACHE_HASH_ARCH}, not one of native, x86-64-v3, x86-64-v4 and none.")
endif()

if(NOT KACHE_HASH_X86 AND NOT KACHE_HASH_ARCH STREQUAL "none")
    if(KACHE_HASH_TOP_LEVEL)
        set(KACHE_HASH_ARCH_CHECK FATAL_ERROR)
    else()
        set(KACHE_HASH_ARCH_CHECK WARNING)
    endif()
    message(${KACHE_HASH_ARCH_CHECK} "kache-hash's compiled targets need an x86-64 CPU with AVX2 for now, and this build "
            "targets ${CMAKE_SYSTEM_PROCESSOR}; with KACHE_HASH_ARCH none, they take the toolchain's flags alone.")
endif()

# Keep branches from crossing or ending on 32-byte boundaries: Skylake-derived cores (such as Cascade Lake) can't serve
# those from the decoded-uop cache under their JCC-erratum microcode. Other cores pay only for the padding. On Cascade
# Lake it gains 1-4% in builds and up to 9% in negative queries.
set(JCC_PAD_FLAGS "")
if(KACHE_HASH_X86 AND NOT KACHE_HASH_ARCH STREQUAL "none")
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(JCC_PAD_FLAG -mbranches-within-32B-boundaries)
    else()
        set(JCC_PAD_FLAG -Wa,-mbranches-within-32B-boundaries)
    endif()

    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag(${JCC_PAD_FLAG} KACHE_HASH_HAS_JCC_PAD)
    if(KACHE_HASH_HAS_JCC_PAD)
        set(JCC_PAD_FLAGS ${JCC_PAD_FLAG})
    endif()
endif()

# Warning flags: https://gcc.gnu.org/onlinedocs/gcc/Warning-Options.html.
set(WARNING_FLAGS -Wall -Wextra -Wpedantic)

# Debug flags: https://gcc.gnu.org/onlinedocs/gcc/Debugging-Options.html.
set(DEBUG_FLAGS -ggdb)

# Address sanitizer flags: https://github.com/google/sanitizers/wiki/AddressSanitizer.
set(ASAN_FLAGS -fsanitize=address -fno-omit-frame-pointer)

# Undefined behavior sanitizer flags: https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html. Some are Clang's
# only.
set(UBSAN_FLAGS -fsanitize=undefined
                -fsanitize=float-divide-by-zero #-fsanitize=unsigned-integer-overflow
                -fsanitize=implicit-conversion -fsanitize=local-bounds -fsanitize=nullability)
if(UBSAN AND NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(FATAL_ERROR "UBSAN takes Clang: configure with -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++.")
endif()

# Compile- and link-flags of the compiled targets.
set(COMPILE_FLAGS ${WARNING_FLAGS} -funroll-loops ${ARCH_FLAGS} ${JCC_PAD_FLAGS})
set(LINK_FLAGS "")

if(DEBUG_INFO)
    list(APPEND COMPILE_FLAGS ${DEBUG_FLAGS})
endif()

if(ASAN)
    list(APPEND COMPILE_FLAGS ${ASAN_FLAGS} ${DEBUG_FLAGS})
    list(APPEND LINK_FLAGS ${ASAN_FLAGS})
endif()

if(UBSAN)
    list(APPEND COMPILE_FLAGS ${UBSAN_FLAGS} ${DEBUG_FLAGS})
    list(APPEND LINK_FLAGS ${UBSAN_FLAGS})
endif()


if(KACHE_HASH_TOP_LEVEL)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)   # Better `clangd` scan.
endif()


# Required packages / libraries from the system.
find_package(Threads REQUIRED)

# Installation directories, which the targets' interfaces name.
include(GNUInstallDirs)


# The table as a header-only library, for other projects to link. It carries what its headers need: C++20, threads, and
# the include path; the headers set their own definitions. Optimization and ISA flags (AVX2 at least) are the linking
# project's own.
add_library(kache_hash INTERFACE)
add_library(kache_hash::kache_hash ALIAS kache_hash)
target_include_directories(kache_hash INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                                                $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_compile_features(kache_hash INTERFACE cxx_std_20)
target_link_libraries(kache_hash INTERFACE Threads::Threads)

# Adds the compiled library `target` of the sources `ARGN`, with its alias `kache_hash::<target>`: static, or shared
# with `BUILD_SHARED_LIBS`, exporting the symbols of its interface only, and position-independent either way. Its
# sources take kache-hash's flags, and Release's in a build without a type, as a parent project's may be; its interface
# takes neither, nor C++20.
function(kache_hash_add_library target)
    add_library(${target} ${ARGN})
    add_library(kache_hash::${target} ALIAS ${target})
    target_include_directories(${target} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                                                $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
    target_compile_options(${target} PRIVATE ${COMPILE_FLAGS} $<$<STREQUAL:$<CONFIG>,>:-O3>)
    target_compile_definitions(${target} PRIVATE $<$<STREQUAL:$<CONFIG>,>:NDEBUG>)
    target_link_options(${target} PRIVATE ${LINK_FLAGS})
    target_link_libraries(${target} PRIVATE kache_hash PUBLIC Threads::Threads)
    set_target_properties(${target} PROPERTIES CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON
                                               POSITION_INDEPENDENT_CODE ON VERSION ${PROJECT_VERSION}
                                               SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
endfunction()

# The k-mer counting library, in `counter/`.
if(KACHE_HASH_BUILD_COUNTER)
    add_subdirectory(counter)
endif()

# The C interface library, in `capi/`.
if(KACHE_HASH_BUILD_CAPI)
    add_subdirectory(capi)
endif()

# The tests, in `test/`, for CTest: in top-level builds only.
if(KACHE_HASH_BUILD_TESTS AND KACHE_HASH_TOP_LEVEL)
    enable_testing()
    add_subdirectory(test)
elseif(KACHE_HASH_BUILD_TESTS)
    message(WARNING "KACHE_HASH_BUILD_TESTS takes effect in top-level builds only.")
endif()


# Installation: the headers, the libraries, the CMake package `kache_hash` with its components `capi` and `counter`, and
# the pkg-config files `kache_hash` and `kache_capi`.
if(KACHE_HASH_INSTALL)
    include(CMakePackageConfigHelpers)
    set(KACHE_HASH_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/kache_hash)

    set(KACHE_HASH_HEADER_EXCLUDES "")
    if(NOT TARGET kache_capi)
        list(APPEND KACHE_HASH_HEADER_EXCLUDES PATTERN kache_capi.h EXCLUDE)
    endif()
    if(NOT TARGET kache_counter)
        list(APPEND KACHE_HASH_HEADER_EXCLUDES PATTERN Kmer_Set_Counter.hpp EXCLUDE PATTERN Fastx_Reader.hpp EXCLUDE)
    endif()
    install(DIRECTORY include/kache_hash DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} ${KACHE_HASH_HEADER_EXCLUDES})
    install(FILES LICENSE DESTINATION ${CMAKE_INSTALL_DOCDIR})
    install(FILES include/kache_hash/xxHash/LICENSE DESTINATION ${CMAKE_INSTALL_DOCDIR} RENAME LICENSE.xxHash)

    install(TARGETS kache_hash EXPORT kache_hash_targets)
    install(EXPORT kache_hash_targets NAMESPACE kache_hash:: DESTINATION ${KACHE_HASH_CMAKE_DIR}
            FILE kache_hashTargets.cmake)
    foreach(component capi counter)
        if(TARGET kache_${component})
            install(TARGETS kache_${component} EXPORT kache_hash_${component}_targets)
            install(EXPORT kache_hash_${component}_targets NAMESPACE kache_hash:: DESTINATION ${KACHE_HASH_CMAKE_DIR}
                    FILE kache_hash_${component}Targets.cmake)
        endif()
    endforeach()

    configure_package_config_file(cmake/kache_hashConfig.cmake.in ${PROJECT_BINARY_DIR}/kache_hashConfig.cmake
                                  INSTALL_DESTINATION ${KACHE_HASH_CMAKE_DIR})
    write_basic_package_version_file(${PROJECT_BINARY_DIR}/kache_hashConfigVersion.cmake
                                     COMPATIBILITY SameMinorVersion)    # A 0.x release may change the interfaces.
    install(FILES ${PROJECT_BINARY_DIR}/kache_hashConfig.cmake ${PROJECT_BINARY_DIR}/kache_hashConfigVersion.cmake
            DESTINATION ${KACHE_HASH_CMAKE_DIR})

    # The pkg-config files are relocatable: their prefix is relative to their own directory.
    foreach(dir INCLUDEDIR LIBDIR)
        if(IS_ABSOLUTE "${CMAKE_INSTALL_${dir}}")
            set(KACHE_HASH_PC_${dir} "${CMAKE_INSTALL_${dir}}")
        else()
            set(KACHE_HASH_PC_${dir} "\${prefix}/${CMAKE_INSTALL_${dir}}")
        endif()
    endforeach()

    file(RELATIVE_PATH KACHE_HASH_PC_PREFIX ${CMAKE_INSTALL_FULL_DATADIR}/pkgconfig ${CMAKE_INSTALL_PREFIX})
    string(REGEX REPLACE "/$" "" KACHE_HASH_PC_PREFIX "${KACHE_HASH_PC_PREFIX}")
    configure_file(cmake/kache_hash.pc.in ${PROJECT_BINARY_DIR}/kache_hash.pc @ONLY)
    install(FILES ${PROJECT_BINARY_DIR}/kache_hash.pc DESTINATION ${CMAKE_INSTALL_DATADIR}/pkgconfig)

    if(TARGET kache_capi)
        # A static library's dependencies, the C++ runtime and threads, are its consumers' to link too.
        if(APPLE)
            set(KACHE_CAPI_PC_DEPS "-lc++")
        else()
            set(KACHE_CAPI_PC_DEPS "-lstdc++ -lm -pthread")
        endif()
        if(BUILD_SHARED_LIBS)
            set(KACHE_CAPI_PC_LIBS "")
            set(KACHE_CAPI_PC_LIBS_PRIVATE "${KACHE_CAPI_PC_DEPS}")
        else()
            set(KACHE_CAPI_PC_LIBS " ${KACHE_CAPI_PC_DEPS}")
            set(KACHE_CAPI_PC_LIBS_PRIVATE "")
        endif()

        file(RELATIVE_PATH KACHE_HASH_PC_PREFIX ${CMAKE_INSTALL_FULL_LIBDIR}/pkgconfig ${CMAKE_INSTALL_PREFIX})
        string(REGEX REPLACE "/$" "" KACHE_HASH_PC_PREFIX "${KACHE_HASH_PC_PREFIX}")
        configure_file(cmake/kache_capi.pc.in ${PROJECT_BINARY_DIR}/kache_capi.pc @ONLY)
        install(FILES ${PROJECT_BINARY_DIR}/kache_capi.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
    endif()
endif()


# The rest is for the benchmarks and the tools, in top-level builds only: in a subproject, they would set the parent's
# build type, write into the parent's source tree, and fetch the baselines.
if(NOT (KACHE_HASH_TOP_LEVEL AND KACHE_HASH_BUILD_BENCHMARKS))
    if(KACHE_HASH_TOP_LEVEL AND (CBL OR SSHASH OR LSWEEP))
        message(WARNING "CBL, SSHASH and LSWEEP take effect with KACHE_HASH_BUILD_BENCHMARKS on.")
    elseif(NOT KACHE_HASH_TOP_LEVEL AND KACHE_HASH_BUILD_BENCHMARKS)
        message(WARNING "KACHE_HASH_BUILD_BENCHMARKS takes effect in top-level builds only.")
    endif()
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


add_subdirectory(benchmark)
add_subdirectory(tools)
