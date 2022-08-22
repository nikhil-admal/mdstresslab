cmake_minimum_required(VERSION 3.5)

enable_testing()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(grid_cauchy.data grid_cauchy.data)
configure_file(grid_pk1.data grid_pk1.data)
configure_file(config.data config.data)
add_executable(testMls testMls.cpp)
target_link_libraries(testMls MDStressLab++)

add_test(TestMls testMls)
