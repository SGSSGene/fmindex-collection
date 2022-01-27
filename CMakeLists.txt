cmake_minimum_required (VERSION 3.8)


# check if this is the main project
if (NOT DEFINED PROJECT_NAME)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts -Wall -Werror")
endif ()

project(fmindex-collection)

find_package(PkgConfig)

pkg_check_modules(DIVSUFSORT64 REQUIRED libdivsufsort64)

#add_subdirectory(lib/sdsl-lite)
add_subdirectory(lib/Catch2)
add_subdirectory(lib/cereal)

add_subdirectory(src/fmindex-collection)
add_subdirectory(src/test_fmindex-collection)
add_subdirectory(src/search_schemes)
add_subdirectory(src/test_search_schemes)

# example executable
add_executable(example
    src/example/main.cpp
    src/example/utils/utils.cpp
)

target_link_libraries(example
    PRIVATE
    fmindex-collection
    search_schemes
    fmt
    procps
)

# easyExample executable
add_executable(easyExample
    src/easyExample/main.cpp
)
target_link_libraries(easyExample
    PRIVATE
    fmindex-collection
    search_schemes
    fmt
    procps
    divsufsort64
)

