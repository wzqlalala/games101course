cmake_minimum_required(VERSION 3.10)
project(Rasterizer)

set(OpenCV_DIR "D:/Program Files/opencv/build/x64/vc15/lib")
find_package(OpenCV REQUIRED)

set(CMAKE_CXX_STANDARD 17)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)
include_directories(/usr/local/include)

add_executable(Rasterizer main.cpp rasterizer.hpp rasterizer.cpp global.hpp Triangle.hpp Triangle.cpp)
target_link_libraries(Rasterizer ${OpenCV_LIBRARIES})
