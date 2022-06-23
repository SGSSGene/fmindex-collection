cmake_minimum_required (VERSION 3.8)

project(test_libsaispp)

add_executable(${PROJECT_NAME}
    sais1.cpp
    sais2.cpp
)
target_link_libraries(${PROJECT_NAME}
    libsaispp
)
add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
