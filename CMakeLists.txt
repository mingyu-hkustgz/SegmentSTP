cmake_minimum_required(VERSION 3.23)
project(SegmentSTP)

set(CMAKE_CXX_STANDARD 17)

include_directories(include)

find_package(OpenMP)

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message(FATAL_ERROR "no OpenMP supprot")
endif()

add_definitions(-std=c++17 -O3)
add_subdirectory(src)
add_subdirectory(test)