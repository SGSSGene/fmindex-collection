cmake_minimum_required (VERSION 3.8)

project(search_schemes)

# search scheme generator library
add_library(${PROJECT_NAME} INTERFACE)
add_library(search_schemes::search_schemes ALIAS ${PROJECT_NAME})
target_include_directories(${PROJECT_NAME}
    INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)
