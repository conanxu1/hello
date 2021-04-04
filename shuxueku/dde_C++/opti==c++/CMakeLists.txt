message("Linux")
cmake_minimum_required (VERSION2.6)
project(dde)


aux_source_directory(. SRC_LIST)
add_executable(ddetest ${SRC_LIST})

add_definitions(-std=c++11)
add_definitions(-g -O3 -mavx -lm -lnlopt  -W)

include_directories(/usr/include/eigen3 )





