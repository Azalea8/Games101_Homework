cmake_minimum_required(VERSION 3.10)
project(BezierCurve)

set(OpenCV_DIR D:/Games101_environment/opencv/MinGW-build/install)
find_package(OpenCV REQUIRED)

set(CMAKE_CXX_STANDARD 14)

add_executable(BezierCurve main.cpp)

target_link_libraries(BezierCurve ${OpenCV_LIBRARIES})
