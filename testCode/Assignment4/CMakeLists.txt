cmake_minimum_required(VERSION 3.10)
project(BezierCurve)
set(OpenCV_DIR "D:/Program Files/opencv/build/x64/vc15/lib")
find_package(OpenCV REQUIRED)

set(CMAKE_CXX_STANDARD 14)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)
include_directories(/usr/local/include ./include)
add_executable(BezierCurve main.cpp)

target_link_libraries(BezierCurve ${OpenCV_LIBRARIES})
