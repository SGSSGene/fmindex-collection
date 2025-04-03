# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0
cmake_minimum_required (VERSION 3.25)

project(fmindex-collection)

# fmindex_collection library
add_library(${PROJECT_NAME} INTERFACE)
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_20)
add_library(fmindex-collection::fmindex-collection ALIAS ${PROJECT_NAME})

target_include_directories(${PROJECT_NAME}
    INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(${PROJECT_NAME}
    INTERFACE
    libsais
    cereal::cereal
)

if (FMC_USE_SDSL)
    target_compile_definitions(${PROJECT_NAME} INTERFACE -DFMC_USE_SDSL)
    target_link_libraries(${PROJECT_NAME} INTERFACE sdsl-lite::sdsl-lite)
endif()
