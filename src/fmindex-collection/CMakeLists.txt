cmake_minimum_required (VERSION 3.8)

project(fmindex-collection)

# fmindex_collection library
add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME}
    INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)
target_include_directories(${PROJECT_NAME}
    SYSTEM INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/sdsl-lite/include
)

target_link_libraries(${PROJECT_NAME}
    INTERFACE
)
