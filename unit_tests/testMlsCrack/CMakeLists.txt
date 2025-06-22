enable_testing()

configure_file(grid_ten.data grid_ten.data)
configure_file(config_T.data config_T.data)
add_executable(testMlsCrack testMlsCrack.cpp)
target_link_libraries(testMlsCrack MDStressLab++)

add_test(testMlsCrack testMlsCrack)
