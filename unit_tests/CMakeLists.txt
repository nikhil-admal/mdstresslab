cmake_minimum_required(VERSION 3.5)

enable_testing()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

enable_testing()

add_subdirectory(testGrid)
add_subdirectory(testSpatialHash)
add_subdirectory(testLJ)
add_subdirectory(testSW)


