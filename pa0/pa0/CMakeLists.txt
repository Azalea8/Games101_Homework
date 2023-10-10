cmake_minimum_required (VERSION 2.8.11)
project (Transformation)

find_package(Eigen3 REQUIRED)
include_directories("C:/Program Files (x86)/Eigen3/include/eigen3")

add_executable (Transformation main.cpp)
