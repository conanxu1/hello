message("Linux")
cmake_minimum_required(VERSION 2.6)
project(dde)


aux_source_directory(. SRC_LIST)
add_executable(ddetest ${SRC_LIST})

add_definitions(-std=c++11)


 
message(${USE_EIGEN})



if(${USE_EIGEN})
message("USE_EIGEN")
add_definitions(-DUSE_EIGEN)
include_directories(/usr/include/eigen3 )
endif(${USE_EIGEN})




find_package(OpenCV REQUIRED)
target_link_libraries(ddetest ${OpenCV_LIBS})

add_definitions(-g -O3 -mavx -lm -lnlopt  -W)
