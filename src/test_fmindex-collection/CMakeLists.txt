# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0
cmake_minimum_required (VERSION 3.8)

project(test_fmindex-collection)

# fmindex-collectionunittests
add_executable(${PROJECT_NAME}
    checkBiFMIndex.cpp
    checkDenseBiFMIndex.cpp
    checkDenseReverseFMIndex.cpp
    checkFMIndex.cpp
    checkOccTables.cpp
    checkReverseFMIndex.cpp
    checkReverseIndexSearch.cpp
    checkSearchBacktracking.cpp
    checkSearchPseudo.cpp
    utils.cpp
)
target_link_libraries(${PROJECT_NAME}
    Catch2::Catch2WithMain
    fmindex-collection
)


target_include_directories(${PROJECT_NAME} PRIVATE ${DIVSUFSORT64_INCLUDE_DIRS})
target_compile_options(${PROJECT_NAME} PRIVATE ${DIVSUFSORT64_CFLAGS})
target_link_libraries(${PROJECT_NAME} ${DIVSUFSORT64_LIBRARIES})
add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
