find_package(CURL REQUIRED)

configure_file(output.xyz output.xyz)
# add the executable
add_executable(newtest main.cpp)
target_link_libraries(newtest PRIVATE MDStressLab++ CURL::libcurl)
