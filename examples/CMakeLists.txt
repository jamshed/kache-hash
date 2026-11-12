
# The examples: in C++, counting k-mers, filtering reads by their k-mers in a reference, strand membership, and
# counting reads' k-mers in a set fixed beforehand, through the counter library; and in C, strand membership through
# the C interface. With the tests, CTest runs them on the inputs in `test/data/`, labeled `quick` and `examples`, and
# compares their outputs with those in `expected/`, which `expected/make-expected.py` makes by a plain count of the
# inputs' 31-mers; the Rust crate's examples print the same for theirs.

# Adds the example program `name` of the source `src`, linked to the libraries `ARGN`. It builds with kache-hash's
# flags, and Release's in a build without a type.
function(kache_hash_example name src)
    add_executable(${name} ${src})
    target_compile_options(${name} PRIVATE ${COMPILE_FLAGS} $<$<STREQUAL:$<CONFIG>,>:-O3>)
    target_compile_definitions(${name} PRIVATE $<$<STREQUAL:$<CONFIG>,>:NDEBUG>)
    target_link_options(${name} PRIVATE ${LINK_FLAGS})
    target_link_libraries(${name} PRIVATE ${ARGN})
endfunction()

kache_hash_example(count_kmers count_kmers.cpp kache_hash)
kache_hash_example(filter_reads filter_reads.cpp kache_hash)
kache_hash_example(strand_membership strand_membership.cpp kache_hash)
if(TARGET kache_counter)
    kache_hash_example(count_in_set count_in_set.cpp kache_counter)
endif()
if(TARGET kache_capi)
    kache_hash_example(strand_membership_c strand_membership.c kache_capi)
endif()


if(NOT KACHE_HASH_BUILD_TESTS)    # The tests take top-level builds, as the examples do.
    return()
endif()

set(data ${PROJECT_SOURCE_DIR}/test/data)

# The inputs, decompressed into the build directory.
add_test(NAME examples.inputs WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
         COMMAND sh -c "gzip -dc '${data}/ecoli.fa.gz' > ecoli.fa && gzip -dc '${data}/reads.fq.gz' > reads.fq")
set_tests_properties(examples.inputs PROPERTIES FIXTURES_SETUP examples LABELS "quick;examples")

# Adds the test `examples.<name>`, which runs the shell command `cmd` in the build directory, on 4 threads at most, and
# compares its output with `expected/<expected>`.
function(kache_hash_example_test name expected cmd)
    add_test(NAME examples.${name} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
             COMMAND sh -c "${cmd} > ${name}.out && diff -u '${CMAKE_CURRENT_SOURCE_DIR}/expected/${expected}' ${name}.out")
    set_tests_properties(examples.${name} PROPERTIES FIXTURES_REQUIRED examples LABELS "quick;examples" PROCESSORS 4)
endfunction()

kache_hash_example_test(count_kmers.ecoli count_kmers.ecoli.txt "$<TARGET_FILE:count_kmers> ecoli.fa 4")
kache_hash_example_test(count_kmers.reads count_kmers.reads.txt "$<TARGET_FILE:count_kmers> reads.fq 4")
kache_hash_example_test(count_kmers.edge-cases count_kmers.edge-cases.txt "$<TARGET_FILE:count_kmers> '${data}/edge-cases.fa' 4")
kache_hash_example_test(filter_reads filter_reads.names "$<TARGET_FILE:filter_reads> ecoli.fa reads.fq 0.5 4 | sed -n 's/^>//p'")
kache_hash_example_test(strand_membership strand_membership.tsv "$<TARGET_FILE:strand_membership> ecoli.fa reads.fq 4")
if(TARGET count_in_set)
    kache_hash_example_test(count_in_set count_in_set.txt "$<TARGET_FILE:count_in_set> ecoli.fa reads.fq 4")
endif()
if(TARGET strand_membership_c)
    kache_hash_example_test(strand_membership_c strand_membership.tsv "$<TARGET_FILE:strand_membership_c> ecoli.fa reads.fq")
endif()
