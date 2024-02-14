# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0
cmake_minimum_required (VERSION 3.12)
if (TARGET fmindex-collection::fmindex-collection)
    return()
endif()

CPMAddPackage(
  NAME libsais
  GITHUB_REPOSITORY IlyaGrebnov/libsais
  GIT_TAG v2.7.3
  OPTIONS
    "LIBSAIS_USE_OPENMP ${OpenMP_C_FOUND}"
    "LIBSAIS_BUILD_SHARED_LIB OFF"
  SYSTEM YES
)
if (FMC_USE_SDSL)
    CPMAddPackage(
      NAME sdsl-lite
      GITHUB_REPOSITORY xxsds/sdsl-lite
      GIT_TAG 206a5f725ee54e892d7cf5f17e77aad4cfb31a62
      DOWNLOAD_ONLY TRUE
    )
    if(sdsl-lite_ADDED)
        add_library(sdsl-lite INTERFACE IMPORTED)
        target_include_directories(sdsl-lite SYSTEM INTERFACE "${sdsl-lite_SOURCE_DIR}/include")
        add_library(sdsl-lite::sdsl-lite ALIAS sdsl-lite)
    endif()
endif()

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src/search_schemes ${CMAKE_CURRENT_BINARY_DIR}/search_schemes)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src/fmindex-collection ${CMAKE_CURRENT_BINARY_DIR}/fmindex-collection)
