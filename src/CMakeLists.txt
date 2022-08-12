cmake_minimum_required(VERSION 3.5)

# find the kim api package
find_package(KIM-API 2.2 REQUIRED CONFIG)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(MDStressLab++ SHARED BoxConfiguration.cpp 
                             Configuration.cpp
                             SubConfiguration.cpp
                             neighbor_list.cpp 
                             Exception.cpp 
                             Sphere.cpp
                             Constant.cpp
                             Trigonometric.cpp
                             kim.cpp
                             InteratomicForces.cpp)

target_link_libraries(MDStressLab++ KIM-API::kim-api)
