cmake_minimum_required(VERSION 3.5)

find_package(CURL REQUIRED)
# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(grid_cauchy.data grid_cauchy.data)
configure_file(config.data config.data)
# add the executable
add_executable(crack main.cpp)
target_link_libraries(crack PRIVATE MDStressLab++ CURL::libcurl)
