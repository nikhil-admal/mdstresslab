find_package(CURL REQUIRED)

configure_file(grid_cauchy.data grid_cauchy.data)
configure_file(config.data config.data)
# add the executable
add_executable(crack main.cpp)
target_link_libraries(crack PRIVATE MDStressLab++ CURL::libcurl)
