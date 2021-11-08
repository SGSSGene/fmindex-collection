cmake_minimum_required (VERSION 3.8)

project(search_schemes)

# search scheme generator library
add_library(${PROJECT_NAME}
    generator/h2.cpp
    generator/pigeon.cpp
    generator/suffixFilter.cpp
    generator/zeroOnesZero.cpp
)
target_include_directories(${PROJECT_NAME}
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)
