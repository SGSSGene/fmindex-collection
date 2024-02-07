# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0
cmake_minimum_required (VERSION 3.8)

project(test_fmindex-collection)

# fmindex-collectionunittests
add_executable(${PROJECT_NAME}
    bitvector/checkBitvector.cpp
    checkOccTables.cpp
    checkReverseIndexSearch.cpp
    checkSearchBacktracking.cpp
    checkSearchPseudo.cpp
    fmindex/checkBiFMIndex.cpp
    fmindex/checkBiFMIndexCursor.cpp
    fmindex/checkDenseBiFMIndex.cpp
    fmindex/checkDenseReverseFMIndex.cpp
    fmindex/checkFMIndex.cpp
    fmindex/checkFMIndexCursor.cpp
    fmindex/checkReverseFMIndex.cpp
    fmindex/checkReverseFMIndexCursor.cpp
    rankvector/checkRankVector.cpp
    utils.cpp
)
target_link_libraries(${PROJECT_NAME}
    Catch2::Catch2WithMain
    fmindex-collection
)

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
