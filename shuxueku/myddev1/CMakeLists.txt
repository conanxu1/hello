message("Linux")
project(DDE)


aux_source_directory(. SRC_LIST)
add_executable(ddetest ${SRC_LIST})

# # SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
# # SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)



