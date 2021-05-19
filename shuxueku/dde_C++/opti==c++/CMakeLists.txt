message("Linux")
cmake_minimum_required(VERSION 2.6)
project(dde)


aux_source_directory(. SRC_LIST)
add_executable(ddetest ${SRC_LIST})
add_definitions(-std=c++11)




option(USE_EIGEN "使用eigen库"  ON)
if(USE_EIGEN)
add_definitions(-DUSE_EIGEN)
include_directories(/usr/include/eigen3 )
endif() 
message("\nEigen")
message(${USE_EIGEN})




option(SMALL_SCALE "小规模"  OFF)
if(SMALL_SCALE)
add_definitions(-DSMALL_SCALE)
endif()
message("\nSMALL_SCALE")
message(${SMALL_SCALE})

option(LARGE_SCALE "大规模"  ON)
if(LARGE_SCALE)
add_definitions(-DLARGE_SCALE)
endif()
message("\nLARGE_SCALE")
message(${LARGE_SCALE})





# # # -DLARGE_SCALE=ON

# # # if(${LARGE_SCALE})
# # # message("LARGE_SCALE")
# # # add_definitions(-DLARGE_SCALE)
# # # endif(${LARGE_SCALE})

# # # -DLARGE_SCALE=1


find_package(OpenCV REQUIRED)
target_link_libraries(ddetest ${OpenCV_LIBS})

add_definitions(-g -O3 -mavx -lm -lnlopt  -W)
