# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0
cmake_minimum_required (VERSION 3.12)
if (TARGET fmindex-collection::fmindex-collection)
    return()
endif()

set(LIBSAIS_USE_OPENMP ${OpenMP_C_FOUND})
set(LIBSAIS_BUILD_SHARED_LIB OFF)
CPMAddPackage("gh:IlyaGrebnov/libsais@2.8.4")

if (FMC_USE_SDSL)
    CPMAddPackage("gh:SGSSGene/sdsl-lite@3.0.3-2")
endif()

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src/search_schemes ${CMAKE_CURRENT_BINARY_DIR}/search_schemes)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src/fmindex-collection ${CMAKE_CURRENT_BINARY_DIR}/fmindex-collection)
