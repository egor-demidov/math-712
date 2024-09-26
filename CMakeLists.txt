cmake_minimum_required(VERSION 3.28)
project(math_712)

add_subdirectory(eigen)

set(CMAKE_CXX_STANDARD 17)

add_executable(math_712_hw_2 hw2.cpp)

target_link_libraries(math_712_hw_2 Eigen3::Eigen)
