enable_testing()

configure_file(triclinic.lmp triclinic.lmp)
add_executable(testTriclinicParser testTriclinicParser.cpp)
target_link_libraries(testTriclinicParser MDStressLab++)

add_test(TestTriclinicParser testTriclinicParser)
