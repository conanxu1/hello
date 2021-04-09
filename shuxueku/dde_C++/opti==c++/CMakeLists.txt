message("Linux")
cmake_minimum_required(VERSION 2.6)
project(dde)


aux_source_directory(. SRC_LIST)
add_executable(ddetest ${SRC_LIST})

add_definitions(-std=c++11)


<<<<<<< HEAD



message("${USE_EIGEN}")



set(USE_EIGEN 1)
if($USE_EIGEN)
include_directories(/usr/include/eigen3 )
endif($USE_EIGEN)
=======
#if(USE_EIGEN)
include_directories(/usr/include/eigen3 )
#endif(USE_EIGEN)
>>>>>>> f7cc98bb44b74f5bdada3d85c34dc9ac57e4d5d5



find_package(OpenCV REQUIRED)
target_link_libraries(ddetest ${OpenCV_LIBS})


 
add_definitions(-g -O3 -mavx -lm -lnlopt  -W)
