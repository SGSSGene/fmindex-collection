# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file.
# -----------------------------------------------------------------------------------------------------
cmake_minimum_required (VERSION 3.12)
if (TARGET fmindex-collection::fmindex-collection)
    return()
endif()

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/lib/libsais ${CMAKE_CURRENT_BINARY_DIR}/libsais)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src/search_schemes ${CMAKE_CURRENT_BINARY_DIR}/search_schemes)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src/fmindex-collection ${CMAKE_CURRENT_BINARY_DIR}/fmindex-collection)
