
# Include paths.
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)


# Executables.
add_executable(measure measure.cpp)


# Configs for targets.

target_include_directories(measure PRIVATE ${INCLUDE_DIR})

target_compile_options(measure PRIVATE ${COMPILE_FLAGS})

target_link_libraries(measure PRIVATE ${LINK_FLAGS})


# Target-links.

target_link_libraries(measure PRIVATE jemalloc ${CMAKE_DL_LIBS} ZLIB::ZLIB Threads::Threads)


# Installation path.
install(TARGETS measure
        DESTINATION bin
)
