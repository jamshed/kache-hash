
# Include paths.
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)


# Targets.


# Executables.
add_executable(example examples/example.cpp)

# Benchmark experiment executables.
add_executable(key_distribution experiments/key-distribution.cpp)
add_executable(minimizer_effect experiments/minimizer-effect.cpp)
add_executable(lf_v_key_dist experiments/lf-vs-key-dist.cpp)


# Dependencies

add_dependencies(key_distribution prj_parlaylib)
add_dependencies(minimizer_effect prj_parlaylib)
add_dependencies(lf_v_key_dist prj_parlaylib)


# Configs for targets.

target_include_directories(example PRIVATE ${INCLUDE_DIR})

target_include_directories(key_distribution PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})
target_include_directories(minimizer_effect PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})
target_include_directories(lf_v_key_dist PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})

target_compile_options(example PRIVATE ${COMPILE_FLAGS})

target_compile_options(key_distribution PRIVATE ${COMPILE_FLAGS})
target_compile_options(minimizer_effect PRIVATE ${COMPILE_FLAGS})
target_compile_options(lf_v_key_dist PRIVATE ${COMPILE_FLAGS})

target_link_libraries(example PRIVATE ${LINK_FLAGS})

target_link_libraries(key_distribution PRIVATE ${LINK_FLAGS})
target_link_libraries(minimizer_effect PRIVATE ${LINK_FLAGS})
target_link_libraries(lf_v_key_dist PRIVATE ${LINK_FLAGS})


# Target-links.

target_link_libraries(example PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads)

target_link_libraries(key_distribution PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads)
target_link_libraries(minimizer_effect PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads)
target_link_libraries(lf_v_key_dist PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads)


# Installation path.
install(TARGETS example key_distribution minimizer_effect lf_v_key_dist
        DESTINATION bin
)
