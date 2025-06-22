cmake_minimum_required(VERSION 3.5)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++20")


# set the project name
project(MDStressLab++ VERSION 1.0)

#set(CMAKE_CXX_COMPILER /opt/homebrew/Cellar/llvm/19.1.3/bin/clang++)
#set(CMAKE_C_COMPILER /opt/homebrew/Cellar/llvm/19.1.3/bin/clang)
#set(CMAKE_CXX_COMPILER /usr/local/opt/llvm/bin/clang++)
#set(CMAKE_C_COMPILER /usr/local/opt/llvm/bin/clang)
set(CMAKE_C_COMPILER /opt/homebrew/opt/llvm/bin/clang)
set(CMAKE_CXX_COMPILER /opt/homebrew/opt/llvm/bin/clang++)

OPTION (USE_OpenMP "Use OpenMP to enable <omp.h>" ON)
if (USE_OpenMP)
    set(OpenMP_C_FLAGS "-fopenmp")
    set(OpenMP_CXX_FLAGS "-fopenmp")
    set(OpenMP_C_LIB_NAMES "omp")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY omp)
    find_package(OpenMP REQUIRED)
endif(USE_OpenMP)

if(OPENMP_FOUND)
    #set(OPENMP_INCLUDES "/usr/local/opt/libomp/include")
    #set(OPENMP_LIBRARIES "/usr/local/opt/libomp/lib")
    #set(OPENMP_INCLUDES "/usr/lib/llvm-14/include")
    #set(OPENMP_LIBRARIES "/usr/lib/llvm-14/lib")
    set(OPENMP_INCLUDES "/opt/homebrew/opt/llvm/include")
    set(OPENMP_LIBRARIES "/opt/homebrew/opt/llvm/lib")
    include_directories("${OPENMP_INCLUDES}")
    link_directories("${OPENMP_LIBRARIES}" "${FORTRAN_LIBRARIES}")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OPENMP_FOUND)


# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include/MDStressLab++)

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(unit_tests)
