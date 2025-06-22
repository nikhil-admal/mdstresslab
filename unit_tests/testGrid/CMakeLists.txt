enable_testing()

configure_file(random_read.dat random_read.dat)
add_executable(testGrid testGrid.cpp)
target_link_libraries(testGrid MDStressLab++)
add_test(TestGrid testGrid)
