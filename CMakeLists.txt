# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file.
# -----------------------------------------------------------------------------------------------------
cmake_minimum_required (VERSION 3.16)

# check if this is the main project
if (NOT DEFINED PROJECT_NAME)
    set(CMAKE_CXX_STANDARD 20)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
    if(WIN32)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj /EHsc")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -pedantic")
    endif()
    set(ROOT_PROJECT TRUE)
else()
    set(ROOT_PROJECT FALSE)
endif ()


project(fmindex-collection LANGUAGES CXX
                            DESCRIPTION "fmindex-collection -- Datastructures and Algorithms for (Bi-)FMIndices and Approximate Pattern Matching")

if(NOT WIN32)
    option(FMC_USE_SDSL "Activate occ tables that use SDSL as a backend" ${ROOT_PROJECT})
else()
    option(FMC_USE_SDSL "Activate occ tables that use SDSL as a backend" OFF)
endif()

find_package(OpenMP QUIET)
option(LIBSAIS_USE_OPENMP "" OpenMP_C_FOUND)

include(fmindex_collection-config.cmake)

if (${ROOT_PROJECT})
    find_package (fmt REQUIRED)
    find_package (cereal REQUIRED HINTS lib)

    add_subdirectory(lib/Catch2)

    enable_testing()
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
        fmt::fmt-header-only
        cereal::cereal
    )


    # easyExample executable
    add_executable(easyExample
        src/easyExample/main.cpp
    )
    target_link_libraries(easyExample
        PRIVATE
        fmindex-collection
        search_schemes
        fmt::fmt-header-only
        cereal::cereal
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
