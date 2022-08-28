cmake_minimum_required(VERSION 3.5)

enable_testing()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(grid_cauchy.data grid_cauchy.data)
configure_file(grid_pk1.data grid_pk1.data)
configure_file(config.data config.data)
configure_file(ldad_constant_refReference.stress ldad_constant_refReference.stress)
configure_file(ldad_constant_defReference.stress ldad_constant_defReference.stress)
configure_file(ldad_trigonometric_refReference.stress ldad_trigonometric_refReference.stress)
configure_file(ldad_trigonometric_defReference.stress ldad_trigonometric_defReference.stress)
add_executable(testLDADSW testLDADSW.cpp)
target_link_libraries(testLDADSW MDStressLab++)

add_test(TestLDADSW testLDADSW)
