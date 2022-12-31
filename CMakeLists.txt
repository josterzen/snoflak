cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(SnowflakeGeneration)

set(CMAKE_CXX_STANDARD 11)

set(SRC_FILES
    src/main.cpp
    src/lattice.cpp
)

add_executable(snowflake ${SRC_FILES})
