cmake_minimum_required (VERSION 3.8)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project (fmindex-collection)

#add_subdirectory(lib/sdsl-lite)
add_subdirectory(lib/Catch2)

# search scheme generator library
add_library(search_schemes
    src/search_schemes/expand.cpp
    src/search_schemes/generator/pigeon.cpp
    src/search_schemes/generator/h2.cpp
    src/search_schemes/generator/suffixFilter.cpp
    src/search_schemes/generator/zeroOnesZero.cpp
    src/search_schemes/isValid.cpp
)
target_include_directories(search_schemes
    PUBLIC
    src
)

# search scheme generator unittests
add_executable(test_search_schemes
    src/test_search_schemes/main.cpp
    src/test_search_schemes/isValid.cpp
    src/test_search_schemes/checkGenerators.cpp
)
target_link_libraries(test_search_schemes
    Catch2::Catch2
    search_schemes
)

# fmindex_collection library
add_library(fmindex_collection INTERFACE)
target_include_directories(fmindex_collection
    INTERFACE
    src/fmindex-collection
    lib/sdsl-lite/include
)

# example executable
add_executable(example
    src/example/main.cpp
    src/example/utils/utils.cpp
)
target_link_libraries(example
    PRIVATE
    fmindex_collection
    search_schemes
    fmt)