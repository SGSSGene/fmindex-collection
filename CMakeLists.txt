cmake_minimum_required (VERSION 3.8)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project (fmindex-collection)

#add_subdirectory(lib/sdsl-lite)

# fmindex_collection library
add_library(fmindex_collection
    src/fmindex-collection/oss/expand.cpp
    src/fmindex-collection/oss/generator/pigeon.cpp
    src/fmindex-collection/oss/generator/h2.cpp
    src/fmindex-collection/oss/generator/suffixFilter.cpp
    src/fmindex-collection/oss/generator/zeroOnesZero.cpp
    src/fmindex-collection/oss/isValid.cpp
    src/fmindex-collection/utils/utils.cpp
)
target_include_directories(fmindex_collection
    PUBLIC
    src/fmindex-collection
    lib/sdsl-lite/include
)

# example executable
add_executable(example
    src/example/main.cpp
)
target_link_libraries(example PRIVATE fmindex_collection fmt)
