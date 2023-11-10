# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0
cmake_minimum_required (VERSION 3.12)

if (TARGET cereal::cereal)
    return()
endif()

add_library(cereal INTERFACE)
target_include_directories(cereal INTERFACE SYSTEM
    ${CMAKE_CURRENT_LIST_DIR}/cereal/include
)
add_library(cereal::cereal ALIAS cereal)
