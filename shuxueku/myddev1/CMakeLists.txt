message("Linux")
project(dde)


aux_source_directory(. SRC_LIST)
add_executable(ddetest ${SRC_LIST})


target_link_libraries(dde PRIVATE m)

# # SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
# # SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)



