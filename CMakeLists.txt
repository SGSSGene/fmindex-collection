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

option(FMC_USE_SDSL "Activate occ tables that use SDSL as a backend" ${PROJECT_IS_TOP_LEVEL})

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(FMC_USE_SDSL OFF)
elseif(CMAKE_SIZEOF_VOID_P EQUAL 4)
    set(FMC_USE_SDSL OFF)
endif()

find_package(OpenMP QUIET)
include(cmake/CPM.cmake)
include(fmindex_collection-config.cmake)
if (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(libsais PRIVATE -fp-model=strict -fno-strict-aliasing)
endif()

if (PROJECT_IS_TOP_LEVEL)
    CPMAddPackage("gh:fmtlib/fmt#11.0.1@11.0.1")
    CPMAddPackage("gh:catchorg/Catch@3.6.0")
    CPMAddPackage("gh:SGSSGene/nanobench@4.3.11-2")
    set(BUILD_DOC OFF)
    set(BUILD_SANDBOX OFF)
    set(SKIP_PERFORMANCE_COMPARISON ON)
    CPMAddPackage("gh:USCiLab/cereal@1.3.2")
    CPMAddPackage("gh:boost-ext/reflect@1.2.1")
    if (reflect_ADDED)
        add_library(reflect INTERFACE)
        target_include_directories(reflect SYSTEM INTERFACE ${reflect_SOURCE_DIR})
        add_library(reflect::reflect ALIAS reflect)
    endif()

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
