# Hardware probes, apart from the benchmarks that use them: the measurement of
# the memory roofs of a roofline analysis, and the hardware counters that it and
# the benchmarks count with (`Perf_Group.hpp`).

add_executable(memory_roof memory-roof.cpp)

target_compile_options(memory_roof PRIVATE ${COMPILE_FLAGS})

target_link_libraries(memory_roof PRIVATE ${LINK_FLAGS} Threads::Threads)


# Installation path.
install(TARGETS memory_roof
        DESTINATION bin
)
