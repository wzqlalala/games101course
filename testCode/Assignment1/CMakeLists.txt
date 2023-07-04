cmake_minimum_required(VERSION 3.10)
project(Rasterizer)
set(OpenCV_DIR "D:/Program Files/opencv/build/x64/vc15/lib")

set(CMAKE_CXX_STANDARD 17)

include_directories(/usr/local/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)

add_executable(Rasterizer main.cpp rasterizer.hpp rasterizer.cpp Triangle.hpp Triangle.cpp)
target_link_libraries(Rasterizer ${OpenCV_LIBRARIES})
