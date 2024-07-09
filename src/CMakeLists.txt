
# Include paths.
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)


# Targets.

# Executables.
add_executable(${PROJECT_NAME} main.cpp)


# Configs for targets.
target_include_directories(${PROJECT_NAME} PRIVATE ${INCLUDE_DIR})
target_compile_options(${PROJECT_NAME} PRIVATE ${COMPILE_FLAGS})
target_link_libraries(${PROJECT_NAME} PRIVATE ${LINK_FLAGS})


# Installation path.
install(TARGETS ${PROJECT_NAME}
        DESTINATION bin
        )
