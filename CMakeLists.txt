cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++14")


# set the project name
project(MDStressLab++ VERSION 1.0)

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include/MDStressLab++)

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(unit_tests)
