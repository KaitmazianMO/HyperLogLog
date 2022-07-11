cmake_minimum_required(VERSION 3.16)

project(HLL LANGUAGES C)

add_subdirectory(hll)

add_executable(HyperLogLog main.c)
target_link_libraries(HyperLogLog m)
target_link_libraries(HyperLogLog hll)
