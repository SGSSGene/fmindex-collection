cmake_minimum_required (VERSION 3.8)

# check if this is the main project
if (NOT DEFINED PROJECT_NAME)
    set(CMAKE_CXX_STANDARD 20)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
    if(NOT WIN32)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -pedantic")
    endif()
    set(ROOT_PROJECT TRUE)
else()
    set(ROOT_PROJECT FALSE)
endif ()



project(fmindex-collection)

option(FMC_USE_SDSL "Activate occ tables that use SDSL as a backend" ${ROOT_PROJECT})

add_subdirectory(src/fmindex-collection)
add_subdirectory(src/search_schemes)

add_library(libsais lib/libsais/src/libsais.c lib/libsais/src/libsais16.c lib/libsais/src/libsais64.c)
target_include_directories(libsais PUBLIC lib/libsais/src)

if (${ROOT_PROJECT})
    find_package (fmt REQUIRED)

    add_subdirectory(lib/Catch2)

    enable_testing()
    add_subdirectory(src/test_fmindex-collection)
    add_subdirectory(src/test_search_schemes)

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

    add_executable(run_search_schemes
        src/run_search_schemes/main.cpp
    )

    target_link_libraries(run_search_schemes
        PRIVATE
        search_schemes
        fmt::fmt-header-only
    )

endif()
