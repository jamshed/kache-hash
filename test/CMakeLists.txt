
# The tests, for CTest, under three labels:
# - `quick`: the test programs on the small inputs in `data/` (see `data/make-fixtures.py`), and a C11 check of the C
#   interface's header;
# - `full`: the programs on the datasets of `scripts/fetch-data.sh small` in `KACHE_HASH_TEST_DATA`, if they are there;
# - `consumer`: the consumers of installed and fetched builds of kache-hash, in `consumer/`.
# The programs build with kache-hash's flags, and with their assertions and the library's, whatever the build type.
# The concurrent tests take `KACHE_HASH_TEST_THREADS` threads, which CTest counts against its `-j`.

find_package(ZLIB REQUIRED)

set(KACHE_HASH_TEST_DATA ${PROJECT_SOURCE_DIR}/data CACHE PATH "The datasets of the full tests (scripts/fetch-data.sh small).")
cmake_host_system_information(RESULT KACHE_HASH_CORES QUERY NUMBER_OF_LOGICAL_CORES)
set(KACHE_HASH_TEST_THREADS ${KACHE_HASH_CORES} CACHE STRING "Threads of the concurrent tests.")
set(T ${KACHE_HASH_TEST_THREADS})


# Adds the test program `name` of the sources `ARGN`.
function(kache_hash_test_program name)
    add_executable(${name} ${ARGN})
    target_include_directories(${name} PRIVATE ${PROJECT_SOURCE_DIR}/support)
    target_compile_options(${name} PRIVATE ${COMPILE_FLAGS} -UNDEBUG)
    target_link_options(${name} PRIVATE ${LINK_FLAGS})
    target_link_libraries(${name} PRIVATE kache_hash ZLIB::ZLIB)
endfunction()

kache_hash_test_program(check check.cpp)
kache_hash_test_program(describers describers.cpp)
kache_hash_test_program(oriented oriented.cpp)
kache_hash_test_program(oriented-lanes oriented.cpp)
target_compile_definitions(oriented-lanes PRIVATE KACHE_HASH_LANES=1)  # The lane describer, whatever the ISA.
set(KACHE_HASH_TABLE_TESTS oriented oriented-lanes)
if(TARGET kache_capi)
    kache_hash_test_program(capi capi.cpp)
    target_link_libraries(capi PRIVATE kache_capi)
    list(APPEND KACHE_HASH_TABLE_TESTS capi)
endif()
if(TARGET kache_counter)
    kache_hash_test_program(counter counter.cpp)
    target_link_libraries(counter PRIVATE kache_counter)
endif()


# Adds the test `name`, of the command `ARGN`, labeled `label`, on `processors` CPUs.
function(kache_hash_test name label processors)
    add_test(NAME ${name} COMMAND ${ARGN})
    set_tests_properties(${name} PROPERTIES LABELS ${label} PROCESSORS ${processors})
endfunction()

# Adds the correctness matrix of `check` on the index `index`, with the negatives `negative`, as the tests
# `<prefix>.<kind>.<threads>.cap<cap>.<mode>`, labeled `label`: sets, and maps of 64-, 32- and 8-bit counts;
# sequential, and concurrent at 1 and `T` threads; presized, and resizing from 100,000 and from 10 k-mers; whole
# fragments, raw sequences, and ragged batches; and concurrent at `T` threads with bucket counts rounded up to powers of
# 2, presized and from 10 k-mers, over whole fragments.
function(kache_hash_check_matrix prefix index negative label)
    foreach(cap 0 100000 10)
        foreach(mode stream raw multi)
            foreach(kind set map map32 map8)
                kache_hash_test(${prefix}.seq-${kind}.1.cap${cap}.${mode} ${label} 1
                                check ${index} ${negative} seq-${kind} 1 ${cap} ${mode})
                kache_hash_test(${prefix}.mt-${kind}.1.cap${cap}.${mode} ${label} 1
                                check ${index} ${negative} mt-${kind} 1 ${cap} ${mode})
                kache_hash_test(${prefix}.mt-${kind}.T.cap${cap}.${mode} ${label} ${T}
                                check ${index} ${negative} mt-${kind} ${T} ${cap} ${mode})
            endforeach()
        endforeach()
    endforeach()
    foreach(cap 0 10)
        foreach(kind set map map32 map8)
            kache_hash_test(${prefix}.mt-${kind}.T.cap${cap}.stream.pow2 ${label} ${T}
                            check ${index} ${negative} mt-${kind} ${T} ${cap} stream pow2)
        endforeach()
    endforeach()
endfunction()

# Adds the tests of the oriented operations and the C interface on the index `index`, with the negatives `negative`,
# as `<program>.<prefix>`, labeled `label`.
function(kache_hash_table_tests prefix index negative label)
    foreach(program IN LISTS KACHE_HASH_TABLE_TESTS)
        kache_hash_test(${program}.${prefix} ${label} ${T} ${program} ${index} ${negative} ${T})
    endforeach()
endfunction()


# Quick tests, on the small inputs.
set(fixtures ${CMAKE_CURRENT_SOURCE_DIR}/data)
kache_hash_test(describers quick 1 describers ${fixtures}/ecoli.fa.gz ${fixtures}/chr21.fa.gz ${fixtures}/reads.fq.gz
                ${fixtures}/edge-cases.fa)
if(TARGET kache_counter)
    kache_hash_test(counter quick ${T} counter check ${CMAKE_CURRENT_BINARY_DIR}/counter-data ${T})
endif()
kache_hash_table_tests(ecoli ${fixtures}/ecoli.fa.gz ${fixtures}/chr21.fa.gz quick)
kache_hash_table_tests(chr21 ${fixtures}/chr21.fa.gz ${fixtures}/ecoli.fa.gz quick)
kache_hash_check_matrix(check.ecoli ${fixtures}/ecoli.fa.gz ${fixtures}/chr21.fa.gz quick)
kache_hash_check_matrix(check.chr21 ${fixtures}/chr21.fa.gz ${fixtures}/ecoli.fa.gz quick)
kache_hash_check_matrix(check.edge-cases ${fixtures}/edge-cases.fa ${fixtures}/ecoli.fa.gz quick)
if(TARGET kache_capi)
    kache_hash_test(capi.header-c11 quick 1 ${CMAKE_C_COMPILER} -std=c11 -Wall -Wextra -Wpedantic -Werror -fsyntax-only
                    -x c ${PROJECT_SOURCE_DIR}/include/kache_hash/kache_capi.h)
endif()


# Full tests, on the datasets, if they are there.
set(data ${KACHE_HASH_TEST_DATA})
if(EXISTS ${data}/ecoli.fa.gz AND EXISTS ${data}/chr21.fa.gz AND EXISTS ${data}/SRR5833294.1M.fastq.gz)
    kache_hash_test(full.describers full 1 describers ${data}/ecoli.fa.gz ${data}/chr21.fa.gz
                    ${data}/SRR5833294.1M.fastq.gz)
    kache_hash_table_tests(full.ecoli ${data}/ecoli.fa.gz ${data}/chr21.fa.gz full)
    kache_hash_table_tests(full.chr21 ${data}/chr21.fa.gz ${data}/ecoli.fa.gz full)
    kache_hash_check_matrix(full.check.ecoli ${data}/ecoli.fa.gz ${data}/chr21.fa.gz full)
    kache_hash_check_matrix(full.check.chr21 ${data}/chr21.fa.gz ${data}/ecoli.fa.gz full)
else()
    message(STATUS "The full tests take the datasets of `scripts/fetch-data.sh small` in KACHE_HASH_TEST_DATA (${data}).")
endif()


# The consumers, of the static and the shared installs and of the sources, with this build's compilers.
set(KACHE_HASH_CONSUMER_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER})
if(CMAKE_OSX_ARCHITECTURES)
    list(APPEND KACHE_HASH_CONSUMER_ARGS -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES})
endif()
kache_hash_test(consumers consumer ${KACHE_HASH_CORES} bash ${CMAKE_CURRENT_SOURCE_DIR}/consumer/run.sh ${CMAKE_CURRENT_BINARY_DIR}/consumers
                ${KACHE_HASH_CONSUMER_ARGS})
set_tests_properties(consumers PROPERTIES TIMEOUT 3600)
