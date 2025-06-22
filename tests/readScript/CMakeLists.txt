cmake_minimum_required(VERSION 3.5)

find_package(CURL REQUIRED)
# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add the executable
add_executable(readScript main.cpp)
target_link_libraries(readScript PRIVATE MDStressLab++ CURL::libcurl)
