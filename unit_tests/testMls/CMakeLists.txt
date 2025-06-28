enable_testing()

configure_file(reference.grid reference.grid)
configure_file(mlsPushedReference.stress mlsPushedReference.stress)
configure_file(config.data config.data)
add_executable(testMls testMls.cpp)
target_link_libraries(testMls MDStressLab++)

add_test(TestMls testMls)
