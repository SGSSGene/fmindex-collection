# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0
cmake_minimum_required (VERSION 3.25)


# check if this is the main project
if (NOT DEFINED PROJECT_NAME)
    set(ROOT_PROJECT TRUE)
else ()
    set(ROOT_PROJECT FALSE)
endif ()

project(fmindex-collection LANGUAGES CXX
                            DESCRIPTION "fmindex-collection -- Datastructures and Algorithms for (Bi-)FMIndices and Approximate Pattern Matching")

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    option(FMC_USE_SDSL "Activate occ tables that use SDSL as a backend" OFF)
else()
    option(FMC_USE_SDSL "Activate occ tables that use SDSL as a backend" ${PROJECT_IS_TOP_LEVEL})
endif()

find_package(OpenMP QUIET)
include(cmake/CPM.cmake)
CPMAddPackage("gh:SGSSGene/WrappedCPMAddPackage@0.2.7")
include(fmindex_collection-config.cmake)
if (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(libsais PRIVATE -fp-model=strict -fno-strict-aliasing)
endif()

if (PROJECT_IS_TOP_LEVEL)
    CPMAddPackage("gh:fmtlib/fmt#10.1.1@10.1.1")
    CPMAddPackage("gh:catchorg/Catch@3.4.0")
    CPMAddPackage("gh:martinus/nanobench@4.3.11")
    WrappedCPMAddPackage("gh:USCiLab/cereal@1.3.2")

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
        fmindex-collection::fmindex-collection
        fmt::fmt-header-only
        cereal::cereal
    )


    # easyExample executable
    add_executable(easyExample
        src/easyExample/main.cpp
    )
    target_link_libraries(easyExample
        PRIVATE
        fmindex-collection::fmindex-collection
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

    add_subdirectory(src/fmindex-collection-stats)
    add_subdirectory(src/test_header)

    if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        target_compile_options(fmindex-collection INTERFACE /bigobj /EHsc -DNOMINMAX)
    else()
        target_compile_options(fmindex-collection INTERFACE -Wall -Werror -Wpedantic -Wextra -fPIC)
    endif()
endif()
