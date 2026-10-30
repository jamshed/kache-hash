
add_library(kache_counter STATIC Fastx_Reader.cpp)
add_library(kache_hash::kache_counter ALIAS kache_counter)

target_include_directories(kache_counter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(kache_counter PRIVATE ${COMPILE_FLAGS} $<$<STREQUAL:$<CONFIG>,>:-O3>)
target_compile_definitions(kache_counter PRIVATE $<$<STREQUAL:$<CONFIG>,>:NDEBUG>)
target_link_libraries(kache_counter PRIVATE kache_hash PUBLIC Threads::Threads)
