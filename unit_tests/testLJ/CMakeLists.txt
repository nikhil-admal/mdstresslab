cmake_minimum_required(VERSION 3.5)

enable_testing()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(grid_cauchy.data grid_cauchy.data)
configure_file(config.data config.data)
configure_file(hardy1Reference.stress hardy1Reference.stress)
configure_file(hardy2Reference.stress hardy2Reference.stress)
configure_file(hardy3Reference.stress hardy3Reference.stress)
configure_file(hardy4Reference.stress hardy4Reference.stress)
configure_file(hardyRandomCauchyReference.stress hardyRandomCauchyReference.stress)
configure_file(hardyRandomPiolaReference.stress hardyRandomPiolaReference.stress)
add_executable(testLJ testLJ.cpp)
target_link_libraries(testLJ MDStressLab++)

add_test(TestLJ testLJ)
