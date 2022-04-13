set(CMAKE_C_COMPILER "/usr/local/bin/gcc-11")
set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-11")

project (HODLRlib)
cmake_minimum_required (VERSION 3.12)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(DTYPE_FLAG "-DUSE_DOUBLE")#DUSE_COMPLEX64
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DTYPE_FLAG} --std=c++17 -O3 -g -fopenmp")
set(EIGEN_PATH "/usr/local/include/eigen3")
set(HOME_PATH "~/Documents/GitHub/Hmatrix")


################ ADDING OF THE NECESSARY HEADER FOLDERS ####################
# Adding folders which contain necessary header files:
include_directories("${EIGEN_PATH}"
                    "${PROJECT_SOURCE_DIR}/include"
                    "${HOME_PATH}/include"
                   )


################ MENTIONING C++ FILES TO LINK WITH HEADERS + FILE TO BUILD ####################

add_library(HODLRlib
    ${PROJECT_SOURCE_DIR}/src/HODLR_Matrix.cpp
    ${PROJECT_SOURCE_DIR}/src/LowRank.cpp
    ${PROJECT_SOURCE_DIR}/src/HODLR_Node.cpp
    ${PROJECT_SOURCE_DIR}/src/HODLR.cpp
    ${PROJECT_SOURCE_DIR}/src/HODLR_NonSPD.cpp
    ${PROJECT_SOURCE_DIR}/src/HODLR_SPD.cpp
    ${PROJECT_SOURCE_DIR}/src/KDTree.cpp
   )
