cmake_minimum_required(VERSION 3.10)
project(Rasterizer)

set(OpenCV_DIR D:/Games101_environment/opencv/MinGW-build/install)
find_package(OpenCV REQUIRED)

set(CMAKE_CXX_STANDARD 17)

include_directories("C:/Program Files (x86)/Eigen3/include/eigen3")

add_executable(Rasterizer main.cpp rasterizer.hpp rasterizer.cpp Triangle.hpp Triangle.cpp)
target_link_libraries(Rasterizer ${OpenCV_LIBRARIES})
