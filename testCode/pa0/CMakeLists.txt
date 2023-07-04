cmake_minimum_required (VERSION 3.24.2)
project (Transformation)

list(APPEND CMAKE_PREFIX_PATH "../") # 安装路径前缀
include_directories(EIGEN3_INCLUDE_DIR)

add_executable (Transformation main.cpp)
