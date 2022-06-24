cmake_minimum_required (VERSION 3.8)

# check if this is the main project
if (NOT DEFINED PROJECT_NAME)
    set(CMAKE_CXX_STANDARD 20)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -pedantic")
    set(ROOT_PROJECT TRUE)
else()
    set(ROOT_PROJECT FALSE)
endif ()



project(fmindex-collection)

option(FMC_USE_SDSL "Activate occ tables that use SDSL as a backend" ${ROOT_PROJECT})

add_subdirectory(src/fmindex-collection)
add_subdirectory(src/search_schemes)
add_subdirectory(src/libsaispp)

if (${ROOT_PROJECT})
    find_package (fmt REQUIRED)

    add_subdirectory(lib/Catch2)

    enable_testing()
    add_subdirectory(src/test_fmindex-collection)
    add_subdirectory(src/test_search_schemes)
    add_subdirectory(src/test_libsaispp)

    # example executable
    add_executable(example
        src/example/main.cpp
        src/example/utils/utils.cpp
    )

    target_include_directories(example
        SYSTEM PRIVATE
        "${CMAKE_CURRENT_LIST_DIR}/lib/cereal/include"
    )
    target_link_libraries(example
        PRIVATE
        fmindex-collection
        search_schemes
        fmt::fmt-header-only
    )


    # easyExample executable
    add_executable(easyExample
        src/easyExample/main.cpp
    )
    target_include_directories(easyExample
        SYSTEM PRIVATE
        "${CMAKE_CURRENT_LIST_DIR}/lib/cereal/include"
    )
    target_link_libraries(easyExample
        PRIVATE
        fmindex-collection
        search_schemes
        fmt::fmt-header-only
    )
endif()
