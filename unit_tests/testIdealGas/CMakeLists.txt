enable_testing()

configure_file(idealGas.lmp idealGas.lmp)
configure_file(in.ideal_gas in.ideal_gas)
add_executable(testIdealGas testIdealGas.cpp)
target_link_libraries(testIdealGas MDStressLab++)

add_test(TestIdealGas testIdealGas)
