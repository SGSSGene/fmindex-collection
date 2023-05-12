# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file.
# -----------------------------------------------------------------------------------------------------
cmake_minimum_required (VERSION 3.12)

if (TARGET cereal::cereal)
    return()
endif()

add_library(cereal INTERFACE)
target_include_directories(cereal INTERFACE SYSTEM
    ${CMAKE_CURRENT_LIST_DIR}/cereal/include
)
add_library(cereal::cereal ALIAS cereal)
