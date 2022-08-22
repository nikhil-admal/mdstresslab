cmake_minimum_required(VERSION 3.5)

enable_testing()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(grid_ten.data grid_ten.data)
configure_file(config_T.data config_T.data)
add_executable(testMlsCrack testMlsCrack.cpp)
target_link_libraries(testMlsCrack MDStressLab++)

add_test(testMlsCrack testMlsCrack)
