cmake_minimum_required (VERSION 3.8)


# check if this is the main project
if (NOT DEFINED PROJECT_NAME)
    set(CMAKE_CXX_STANDARD 20)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts -Wall -Werror")
    set(ROOT_PROJECT TRUE)
endif ()

project(fmindex-collection)


add_subdirectory(src/fmindex-collection)
add_subdirectory(src/search_schemes)

if (DEFINED ROOT_PROJECT)
    add_subdirectory(lib/Catch2)
    add_subdirectory(lib/cereal)

    add_subdirectory(src/test_fmindex-collection)
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
        cereal
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
        cereal
        fmt
        procps
    )
endif()
