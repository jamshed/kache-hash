# Hardware probes, apart from the benchmarks that use them: the measurement of
# the memory roofs of a roofline analysis, and the hardware counters that it and
# the benchmarks count with (`Perf_Group.hpp`).

target_compile_options(memory_roof PRIVATE ${COMPILE_FLAGS})

target_link_libraries(memory_roof PRIVATE ${LINK_FLAGS})
