cmake_minimum_required(VERSION 3.5)

# find the kim api package
find_package(KIM-API 2.0 REQUIRED)
add_kim_api_imported_library()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(MDStressLab++ SHARED BoxConfiguration.cpp 
                             Configuration.cpp
                             SubConfiguration.cpp
                             neighbor_list.cpp 
                             Exception.cpp 
                             Hardy.cpp
                             kim.cpp
                             InteratomicForces.cpp)

target_link_libraries(MDStressLab++ kim-api)
