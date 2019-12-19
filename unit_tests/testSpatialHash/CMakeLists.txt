cmake_minimum_required(VERSION 3.5)

enable_testing()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(testSpatialHash testSpatialHash.cpp)
target_link_libraries(testSpatialHash MDStressLab++)
add_test(TestSpatialHash  testSpatialHash)
