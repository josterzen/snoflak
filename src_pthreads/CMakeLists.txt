cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(SnowflakeGeneration)

set(CMAKE_CXX_STANDARD 11)

set(SRC_FILES
    main.cpp
    lattice.cpp
)

add_executable(snowflake ${SRC_FILES})

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(snowflake PRIVATE Threads::Threads)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(snowflake PUBLIC OpenMP::OpenMP_CXX)
endif()
