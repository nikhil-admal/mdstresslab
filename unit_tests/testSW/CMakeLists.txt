enable_testing()

configure_file(grid_cauchy.data grid_cauchy.data)
configure_file(config.data config.data)
configure_file(hardy1Reference.stress hardy1Reference.stress)
configure_file(hardy2Reference.stress hardy2Reference.stress)
configure_file(hardy3Reference.stress hardy3Reference.stress)
configure_file(hardy4Reference.stress hardy4Reference.stress)
configure_file(hardyRandomCauchyReference.stress hardyRandomCauchyReference.stress)
configure_file(hardyRandomPiolaReference.stress hardyRandomPiolaReference.stress)
add_executable(testSW testSW.cpp)
target_link_libraries(testSW MDStressLab++)

add_test(TestSW testSW)
