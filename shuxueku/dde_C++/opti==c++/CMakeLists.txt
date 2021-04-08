message("Linux")
cmake_minimum_required(VERSION 2.6)
project(dde)


aux_source_directory(. SRC_LIST)
add_executable(ddetest ${SRC_LIST})

add_definitions(-std=c++11)

include_directories(/usr/include/eigen3 )




find_package(OpenCV REQUIRED)
target_link_libraries(ddetest ${OpenCV_LIBS})


#set(CROSS_COMPILE 1)

if(CROSS_COMPILE)
set(CMAKE_C_COMPILER "arm-linux-gnueabi-gcc")
set(CMAKE_CXX_COMPILER "arm-linux-gnueabi-g++")
add_definitions(-g -O3 -marm -lm -lnlopt  -W)
endif(CROSS_COMPILE)


set(ORD 1)
if(ORD)
add_definitions(-g -O3 -mavx -lm -lnlopt  -W)
endif(ORD)
