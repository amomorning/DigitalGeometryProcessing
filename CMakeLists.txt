cmake_minimum_required (VERSION 3.8)
project("DigitalGeometryProcessing")
set(CMAKE_BUILD_TYPE "Release")
set(LIB_DIR "C:/Users/amomorning/Documents/external/x64")

# Eigen
set(EIGEN_DIR "${LIB_DIR}/eigen3")
include_directories("${EIGEN_DIR}")

# SurfaceMesh
set(SURFACEMESH_DIR "${LIB_DIR}/surface_mesh")
include_directories("${SURFACEMESH_DIR}/include")
link_libraries("${SURFACEMESH_DIR}/lib/surface_mesh.lib")


include_directories("include")
# src - deformation
add_executable(Deformation "Deformation/main.cpp")

# src - Mapping
add_executable(Mapping "Mapping/main.cpp" "Mapping/laplace.cpp" "Mapping/single_patch_parameterization.cpp" )

#src - test
add_executable(Test "Test/main.cpp")