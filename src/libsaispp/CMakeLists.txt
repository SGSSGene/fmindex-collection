cmake_minimum_required (VERSION 3.8)

project(libsaispp)

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME}
    INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)
target_include_directories(${PROJECT_NAME}
    SYSTEM INTERFACE
)

target_link_libraries(${PROJECT_NAME}
    INTERFACE
)
