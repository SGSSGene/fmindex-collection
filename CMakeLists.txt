# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0
cmake_minimum_required (VERSION 3.25)

project(fmindex-collection LANGUAGES CXX C
                            DESCRIPTION "fmindex-collection -- Datastructures and Algorithms for (Bi-)FMIndices and Approximate Pattern Matching")

option(FMC_USE_SDSL  "Include SDSL based strings with rank support" ${PROJECT_IS_TOP_LEVEL})
option(FMC_USE_PASTA "Include PASTA (required for bitvector benchmarks)" ${PROJECT_IS_TOP_LEVEL})
option(FMC_USE_SUX   "Include SUX (required for bitvector benchmarks)" ${PROJECT_IS_TOP_LEVEL})
option(FMC_USE_AWFMINDEX  "Include AWFMIndex (only used in string benchmark)" ${PROJECT_IS_TOP_LEVEL})

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(FMC_USE_SDSL OFF)
    set(FMC_USE_PASTA OFF)
    set(FMC_USE_SUX OFF)
    if (PROJECT_IS_TOP_LEVEL)
        add_compile_options("/EHsc")
    endif()
elseif (EMSCRIPTEN)
    set(FMC_USE_SDSL OFF)
    set(FMC_USE_PASTA OFF)
    set(FMC_USE_SUX OFF)
    set(FMC_USE_AWFMINDEX OFF)
endif()

find_package(OpenMP QUIET)
include(cmake/CPM.cmake)
CPMAddPackage("gh:SGSSGene/cpm.dependencies@1.0.0")
CPMLoadDependenciesFile("${CMAKE_CURRENT_SOURCE_DIR}/cpm.dependencies")
include(fmindex_collection-config.cmake)
if (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(libsais PRIVATE -fp-model=strict -fno-strict-aliasing)
endif()

if (PROJECT_IS_TOP_LEVEL)
    enable_testing()
    add_subdirectory(src/test_fmindex-collection)

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

    # example executable
    add_executable(search_scheme_generator
        src/search_scheme_generator/main.cpp
    )

    target_link_libraries(search_scheme_generator
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
        fmindex-collection::fmindex-collection
        fmt::fmt-header-only
    )

    add_subdirectory(src/fmindex-collection-stats)
    add_subdirectory(src/test_header)
    add_subdirectory(src/docs_examples)

    if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        target_compile_options(fmindex-collection INTERFACE /bigobj -DNOMINMAX)
    else()
        target_compile_options(fmindex-collection INTERFACE -Wall -Werror -Wpedantic -Wextra -fPIC)
    endif()
endif()
