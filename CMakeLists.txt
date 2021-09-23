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
add_library(fmindex_oss
    src/oss/expand.cpp
    src/oss/generator/pigeon.cpp
    src/oss/generator/h2.cpp
    src/oss/generator/suffixFilter.cpp
    src/oss/generator/zeroOnesZero.cpp
    src/oss/isValid.cpp
)
target_include_directories(fmindex_oss
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
    fmindex_oss
    fmt)
