# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file.
# -----------------------------------------------------------------------------------------------------
cmake_minimum_required (VERSION 3.8)

project(fmindex-collection)

# fmindex_collection library
add_library(${PROJECT_NAME} INTERFACE)
add_library(fmindex-collection::fmindex-collection ALIAS ${PROJECT_NAME})

option(FMC_USE_SDSL "Activate occ tables that use SDSL as a backend" ${ROOT_PROJECT})

target_include_directories(${PROJECT_NAME}
    INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(${PROJECT_NAME}
    INTERFACE
    libsais
    search_schemes::search_schemes
)

if (FMC_USE_SDSL)
    target_compile_definitions(${PROJECT_NAME} INTERFACE -DFMC_USE_SDSL)
    target_include_directories(${PROJECT_NAME}
        SYSTEM INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/sdsl-lite/include
    )
endif()
