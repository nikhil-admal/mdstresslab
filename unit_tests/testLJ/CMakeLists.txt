cmake_minimum_required(VERSION 3.5)

enable_testing()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(testLJ testLJ.cpp)
target_link_libraries(testLJ MDStressLab++)

add_test(TestLJ testLJ)
