find_package(CURL REQUIRED)

# add the executable
add_executable(plate main.cpp)
target_link_libraries(plate PRIVATE MDStressLab++ CURL::libcurl)
