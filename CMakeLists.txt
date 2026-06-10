cmake_minimum_required(VERSION 3.5)

# set the project name
project(MDStressLab++ VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


OPTION (USE_OpenMP "Use OpenMP to enable <omp.h>" ON)
if (USE_OpenMP)
    # AppleClang needs explicit hints for OpenMP
    if(APPLE)
        set(OpenMP_CXX_FLAGS "-Xclang -fopenmp")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
    endif()

    find_package(OpenMP COMPONENTS CXX REQUIRED)
endif()


# Include directories
#include_directories(${PROJECT_SOURCE_DIR}/include/MDStressLab++)

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(unit_tests)
