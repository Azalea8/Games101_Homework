cmake_minimum_required(VERSION 3.10)
project(BezierCurve)

set(OpenCV_DIR D:/Code/CLionProjects/Games101_Homework/Games101_environment/OpenCV-MinGW-Build-OpenCV-4.5.2-x64)
find_package(OpenCV REQUIRED)

set(CMAKE_CXX_STANDARD 14)

add_executable(BezierCurve main.cpp)

target_link_libraries(BezierCurve ${OpenCV_LIBRARIES})
