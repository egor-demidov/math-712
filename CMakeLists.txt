cmake_minimum_required(VERSION 3.28)
project(math_712)

add_subdirectory(eigen)#

set(CMAKE_CXX_STANDARD 20)

find_package(OpenMP)

set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Werror -Wuninitialized ${CMAKE_CXX_FLAGS}")
set(CMAKE_C_FLAGS "-O3 -Wall -Wextra -Wuninitialized ${CMAKE_C_FLAGS}")

add_executable(math_712_hw_2_part_1 hw2_part1.cpp)
add_executable(math_712_hw_2_part_2 hw2_part2.cpp)
add_executable(hw2 hw2.c)

if (${OpenMP_FOUND})
    target_link_libraries(math_712_hw_2_part_1 OpenMP::OpenMP_CXX)
    target_link_libraries(math_712_hw_2_part_2 OpenMP::OpenMP_CXX)
    target_link_libraries(hw2 m OpenMP::OpenMP_C)
endif ()

#target_link_libraries(math_712_hw_2_part_1 Eigen3::Eigen)
#target_link_libraries(math_712_hw_2_part_2 Eigen3::Eigen)
