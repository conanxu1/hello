message("Linux")
cmake_minimum_required(VERSION 2.6)
project(dde)


aux_source_directory(. SRC_LIST)
add_executable(ddetest ${SRC_LIST})

add_definitions(-std=c++11)




option(USE_EIGEN "使用eigen库"  OFF)

message(${USE_EIGEN})


if(USE_EIGEN)
add_definitions(-DUSE_EIGEN)
include_directories(/usr/include/eigen3 )
endif()
#-----------------------算值
# # if(${USE_EIGEN})
# # message("USE_EIGEN")
# # add_definitions(-DUSE_EIGEN)
# # include_directories(/usr/include/eigen3 )
# # endif(${USE_EIGEN})


if(${SMALL_SCALE})
message("SMALL_SCALE")
add_definitions(-DSMALL_SCALE)
endif(${SMALL_SCALE})



if(LARGE_SCALE)
message("LARGE_SCALE")
add_definitions(-DLARGE_SCALE)
endif()


# # if(${LARGE_SCALE})
# # message("LARGE_SCALE")
# # add_definitions(-DLARGE_SCALE)
# # endif(${LARGE_SCALE})




find_package(OpenCV REQUIRED)
target_link_libraries(ddetest ${OpenCV_LIBS})

add_definitions(-g -O3 -mavx -lm -lnlopt  -W)
