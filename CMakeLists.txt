cmake_minimum_required(VERSION 2.8.9)
project(example)

add_definitions(-Weffc++)

set(cmake_build_type release)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

file(GLOB SOURCES "*.cpp")
add_executable(example ${SOURCES})
set_property(TARGET example PROPERTY CXX_STANDARD 14)
 
