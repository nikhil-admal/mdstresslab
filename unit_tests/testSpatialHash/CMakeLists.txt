enable_testing()

add_executable(testSpatialHash testSpatialHash.cpp)
target_link_libraries(testSpatialHash MDStressLab++)
add_test(TestSpatialHash  testSpatialHash)
