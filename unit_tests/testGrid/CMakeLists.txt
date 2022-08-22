cmake_minimum_required(VERSION 3.5)

enable_testing()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(random_read.dat random_read.dat)
add_executable(testGrid testGrid.cpp)
target_link_libraries(testGrid MDStressLab++)
add_test(TestGrid testGrid)
