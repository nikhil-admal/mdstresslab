enable_testing()

configure_file(grid_cauchy.data grid_cauchy.data)
configure_file(grid_pk1.data grid_pk1.data)
configure_file(config.data config.data)
configure_file(../testLDADSW/ldad_constant_refReference.stress ldad_constant_refReference.stress)
configure_file(../testLDADSW/ldad_constant_defReference.stress ldad_constant_defReference.stress)
configure_file(../testLDADSW/ldad_trigonometric_refReference.stress ldad_trigonometric_refReference.stress)
configure_file(../testLDADSW/ldad_trigonometric_defReference.stress ldad_trigonometric_defReference.stress)
add_executable(testLDADSWTriclinic testLDADSWTriclinic.cpp)
target_link_libraries(testLDADSWTriclinic MDStressLab++)

add_test(TestLDADSWTriclinic testLDADSWTriclinic)
