cmake_minimum_required(VERSION 3.10)
project(Rasterizer)

set(OpenCV_DIR D:/Games101_environment/OpenCV-MinGW-Build-OpenCV-4.5.2-x64)
find_package(OpenCV REQUIRED)

set(CMAKE_CXX_STANDARD 17)

include_directories("D:/Games101_environment/Eigen3/include/eigen3")

add_executable(Rasterizer main.cpp rasterizer.hpp rasterizer.cpp global.hpp Triangle.hpp Triangle.cpp Texture.hpp Texture.cpp Shader.hpp OBJ_Loader.h)
target_link_libraries(Rasterizer ${OpenCV_LIBRARIES})
#target_compile_options(Rasterizer PUBLIC -Wall -Wextra -pedantic)