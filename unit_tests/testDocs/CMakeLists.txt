enable_testing()

configure_file(grid_cauchy.data grid_cauchy.data)
configure_file(grid_pk1.data grid_pk1.data)
configure_file(config.data config.data)
add_executable(testDocs testDocs.cpp)
target_link_libraries(testDocs MDStressLab++)

add_test(TestDocs testDocs)
