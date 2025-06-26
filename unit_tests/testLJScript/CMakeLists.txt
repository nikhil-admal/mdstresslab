enable_testing()

configure_file(configCurrent.lmp configCurrent.lmp)
configure_file(readScript.in readScript.in)
add_test(TestLJScript ../../src/MDStressLab++_in < readScript.in)
