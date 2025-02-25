# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0
cmake_minimum_required (VERSION 3.8)

project(test_fmindex-collection)

# fmindex-collectionunittests
add_executable(${PROJECT_NAME}
    bitvector/bitset_popcount.cpp
    bitvector/checkBitvector.cpp
    fmindex/checkBiFMIndex.cpp
    fmindex/checkBiFMIndexCursor.cpp
    fmindex/checkCSA.cpp
    fmindex/checkDenseBiFMIndex.cpp
    fmindex/checkDenseReverseFMIndex.cpp
    fmindex/checkFMIndex.cpp
    fmindex/checkFMIndexCursor.cpp
    fmindex/checkLeftBiFMIndexCursor.cpp
    fmindex/checkLeftRBiFMIndexCursor.cpp
    fmindex/checkMerge.cpp
    fmindex/checkRBiFMIndex.cpp
    fmindex/checkRBiFMIndexCursor.cpp
    fmindex/checkReverseFMIndex.cpp
    fmindex/checkReverseFMIndexCursor.cpp
    occtables/checkOccTables.cpp
    rankvector/checkRankVector.cpp
    search/checkReverseIndexSearch.cpp
    search/checkSearchBacktracking.cpp
    search/checkSearchPseudo.cpp
    search/checkLocateFMTree.cpp
    search/checkSearches.cpp
    utils.cpp
    search/checkSearchHammingSM.cpp
)
target_link_libraries(${PROJECT_NAME}
    Catch2::Catch2WithMain
    fmindex-collection
    nanobench::nanobench
    reflect::reflect
    fmt::fmt
)

if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_link_options(${PROJECT_NAME} PUBLIC /STACK:4194304)
endif()

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
