message("Linux")
project(dde)


aux_source_directory(. SRC_LIST)
add_executable(ddetest ${SRC_LIST})

add_definitions(-std=c++11)
add_definitions(-g -O3 -mavx -lm)


include_directories(/usr/include/eigen3 )

# # # set(CMAKE_CXX_STANDARD 11)
# # # target_link_libraries(-lm -lstdc++)

# # SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
# # SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)



