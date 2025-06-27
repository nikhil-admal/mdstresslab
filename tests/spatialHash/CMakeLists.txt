configure_file(NN_full_grid_eps.dat NN_full_grid_eps.dat)
# add the executable
add_executable(spatialHash spatialHash.cpp)
target_link_libraries(spatialHash MDStressLab++)
