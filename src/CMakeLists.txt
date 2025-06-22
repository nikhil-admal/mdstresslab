# find the kim api package
find_package(KIM-API 2.2 REQUIRED CONFIG)

add_library(MDStressLab++ SHARED BoxConfiguration.cpp 
                             Configuration.cpp
                             SubConfiguration.cpp
                             neighbor_list.cpp 
                             Exception.cpp 
                             MethodSphere.cpp
                             Constant.cpp
                             Trigonometric.cpp
                             kim.cpp
                             InteratomicForces.cpp
                             Stencil.cpp
                             Mls.cpp
                             Rigidity.cpp)

target_link_libraries(MDStressLab++ KIM-API::kim-api)

add_executable(MDStressLab++_in readScript.cpp)
target_link_libraries(MDStressLab++_in PRIVATE MDStressLab++)

