cmake_minimum_required(VERSION 3.28)
project(math_712 C)

find_package(OpenMP)

set(CMAKE_C_FLAGS "-O3 -Wall -Wextra -Wuninitialized -Rpass-analysis=loop-vectorize -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -march=native ${CMAKE_C_FLAGS}")

add_executable(hw2 hw2/hw2.c hw2/heat_equation_fwd.c csv_writer.c)

if (${OpenMP_FOUND})
    target_link_libraries(hw2 m OpenMP::OpenMP_C)
else ()
    target_link_libraries(hw2 m)
endif ()

#target_link_libraries(math_712_hw_2_part_1 Eigen3::Eigen)
#target_link_libraries(math_712_hw_2_part_2 Eigen3::Eigen)
