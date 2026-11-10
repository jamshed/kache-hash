# Checks that the version of the headers (`include/kache_hash/version.hpp`), which the CMake project takes, is the
# crate's (`rust/kache-hash-sys/Cargo.toml`) and the citation's (`CITATION.cff`): a release bumps them together.
# Usage: cmake -DROOT=<source-dir> -P cmake/CheckVersion.cmake

file(STRINGS "${ROOT}/include/kache_hash/version.hpp" defs REGEX "^#define KACHE_HASH_VERSION_(MAJOR|MINOR|PATCH) [0-9]+$")
set(headers "")
foreach(part MAJOR MINOR PATCH)
    if(NOT defs MATCHES "KACHE_HASH_VERSION_${part} ([0-9]+)")
        message(FATAL_ERROR "kache_hash/version.hpp lacks KACHE_HASH_VERSION_${part}.")
    endif()
    list(APPEND headers ${CMAKE_MATCH_1})
endforeach()
list(JOIN headers "." headers)

file(STRINGS "${ROOT}/rust/kache-hash-sys/Cargo.toml" crate REGEX "^version = \"[^\"]+\"$" LIMIT_COUNT 1)
string(REGEX REPLACE "^version = \"([^\"]+)\"$" "\\1" crate "${crate}")

file(STRINGS "${ROOT}/CITATION.cff" citation REGEX "^version: .+$" LIMIT_COUNT 1)
string(REGEX REPLACE "^version: (.+)$" "\\1" citation "${citation}")

if(NOT crate STREQUAL headers OR NOT citation STREQUAL headers)
    message(FATAL_ERROR "The versions disagree: the headers' ${headers}, the crate's ${crate}, and CITATION.cff's ${citation}.")
endif()
message(STATUS "The headers, the crate and CITATION.cff are all at ${headers}.")
