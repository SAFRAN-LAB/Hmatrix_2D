set(CMAKE_C_COMPILER "/usr/local/bin/gcc-11")
set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-11")

project(HMAT_test)

cmake_minimum_required (VERSION 3.12)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(EIGEN_PATH "/usr/local/include/eigen3")
set(HOME_PATH "~/Documents/GitHub/Hmatrix")
set(HODLR_PATH "~/Documents/GitHub/Hmatrix/HODLR-master/include")

# add HODLR2D Library
add_subdirectory(HODLR-master)

include_directories("${EIGEN_PATH}"
                    "${HODLR_PATH}"
                    "${HOME_PATH}/include")

set(DTYPE_FLAG "-DCOMPUTE -DUSE_DOUBLE") # -DCOMPUTE -DUSE_COMPLEX64  -DUSE_DOUBLE
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DTYPE_FLAG} --std=c++17 -O3 -g -fopenmp")

# add the executable
add_executable(test_hmat test_ker.cpp)

target_link_libraries(test_hmat PUBLIC HODLRlib)
# add the binary tree to the search path for include files
target_include_directories(test_hmat PUBLIC
                      "${EIGEN_PATH}"
                      "${HODLR_PATH}"
                      "${HOME_PATH}/include")
