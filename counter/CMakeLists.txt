
# The k-mer counting library: a parallel FASTA and FASTQ reader, and a counter
# of reads' k-mers within a set, on the table. Its interface is plain C++, so
# that projects linking it compile none of the table's code and keep their own
# flags. Its own sources take kache-hash's: the warning, optimization and ISA
# flags, and Release's in a build without a type, as a parent project's may be.

add_library(kache_counter STATIC Fastx_Reader.cpp Kmer_Set_Counter.cpp)
add_library(kache_hash::kache_counter ALIAS kache_counter)

target_include_directories(kache_counter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(kache_counter PRIVATE ${COMPILE_FLAGS} $<$<STREQUAL:$<CONFIG>,>:-O3>)
target_compile_definitions(kache_counter PRIVATE $<$<STREQUAL:$<CONFIG>,>:NDEBUG>)
target_link_libraries(kache_counter PRIVATE kache_hash PUBLIC Threads::Threads)
