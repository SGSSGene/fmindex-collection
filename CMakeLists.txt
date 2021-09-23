cmake_minimum_required (VERSION 3.8)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project (fmindex-collection)

#add_subdirectory(lib/sdsl-lite)

# fmindex_collection library
add_library(fmindex_collection INTERFACE)
target_include_directories(fmindex_collection
    INTERFACE
    src/fmindex-collection
    lib/sdsl-lite/include
)

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
