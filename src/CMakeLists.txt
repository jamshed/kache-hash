
# Include paths.
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)


# Targets.

# Libraries.

set(SRC
        Ref_Parser.cpp)

add_library(core STATIC ${SRC})

# Executables.
add_executable(${PROJECT_NAME} main.cpp)


# Configs for targets.

target_include_directories(core PRIVATE ${INCLUDE_DIR} ${EXT_INCLUDE})
target_include_directories(${PROJECT_NAME} PRIVATE ${INCLUDE_DIR})

target_compile_options(core PRIVATE ${COMPILE_FLAGS})
target_compile_options(${PROJECT_NAME} PRIVATE ${COMPILE_FLAGS})

target_link_libraries(core PRIVATE ${LINK_FLAGS})
target_link_libraries(${PROJECT_NAME} PRIVATE ${LINK_FLAGS})


# Target-links.
target_link_libraries(core PRIVATE jemalloc)
# target_link_libraries(core PRIVATE ${CMAKE_DL_LIBS})

target_link_libraries(core PRIVATE ZLIB::ZLIB)

target_link_libraries(core PRIVATE Threads::Threads)

target_link_libraries(${PROJECT_NAME} PRIVATE core)


# Installation path.
install(TARGETS ${PROJECT_NAME}
        DESTINATION bin
        )
