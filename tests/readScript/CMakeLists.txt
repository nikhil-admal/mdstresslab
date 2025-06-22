find_package(CURL REQUIRED)

# add the executable
add_executable(readScript main.cpp)
target_link_libraries(readScript PRIVATE MDStressLab++ CURL::libcurl)
