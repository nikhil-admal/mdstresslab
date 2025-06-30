enable_testing()

configure_file(grid_cauchy.data grid_cauchy.data)
configure_file(config.data config.data)
configure_file(hardyCauchyReference.stress hardyCauchyReference.stress)
configure_file(hardyPiolaReference.stress hardyPiolaReference.stress)
configure_file(virialCauchyReference.stress virialCauchyReference.stress)
configure_file(virialPiolaReference.stress virialPiolaReference.stress)
add_executable(testLJ testLJ.cpp)
target_link_libraries(testLJ MDStressLab++)

add_test(TestLJ testLJ)
