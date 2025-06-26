enable_testing()

configure_file(configCurrent.lmp configCurrent.lmp)
configure_file(testLJScript.in testLJScript.in)
add_test(TestLJScript COMMAND sh -c "../../src/MDStressLab++_in < testLJScript.in")
