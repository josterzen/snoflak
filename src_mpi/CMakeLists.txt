cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(SnowflakeGeneration)

set(CMAKE_CXX_STANDARD 11)

set(SRC_FILES
    main.cpp
    lattice.cpp
)

add_executable(snowflake ${SRC_FILES})

find_package(MPI REQUIRED)
target_link_libraries(snowflake PUBLIC MPI::MPI_CXX)
