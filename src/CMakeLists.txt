
# Include paths.
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)


# Targets.

# Libraries.

set(SRC
        Parser.cpp
        examples/Kmer_Set.cpp
        examples/Minimizer_Scan.cpp
)

add_library(core STATIC ${SRC})

# Executables.
add_executable(bench bench.cpp)
add_executable(main main.cpp)


# Dependencies
add_dependencies(core prj_benchmark)


# Configs for targets.

target_include_directories(core PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})
target_include_directories(bench PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})
target_include_directories(main PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})

target_compile_options(core PRIVATE ${COMPILE_FLAGS})
target_compile_options(bench PRIVATE ${COMPILE_FLAGS})
target_compile_options(main PRIVATE ${COMPILE_FLAGS})

target_link_libraries(core PRIVATE ${LINK_FLAGS})
target_link_libraries(bench PRIVATE ${LINK_FLAGS})
target_link_libraries(main PRIVATE ${LINK_FLAGS})


# Target-links.
target_link_libraries(core PRIVATE jemalloc)
target_link_libraries(core PRIVATE ${CMAKE_DL_LIBS})    # For dynamic shared objects, like `jemalloc`.

target_link_libraries(core PRIVATE ZLIB::ZLIB)

target_link_libraries(core PRIVATE Threads::Threads)

target_link_libraries(core PRIVATE benchmark)

target_link_libraries(core PRIVATE pfm)

target_link_libraries(bench PRIVATE core)
target_link_libraries(bench PRIVATE benchmark)
target_link_libraries(bench PRIVATE nthash)

target_link_libraries(main PRIVATE core)
target_link_libraries(main PRIVATE benchmark)
target_link_libraries(main PRIVATE nthash)


# Installation path.
install(TARGETS bench main
        DESTINATION bin
        )
