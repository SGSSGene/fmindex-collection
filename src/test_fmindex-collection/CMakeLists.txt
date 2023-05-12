# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file.
# -----------------------------------------------------------------------------------------------------
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
    main.cpp
)
target_link_libraries(${PROJECT_NAME}
    Catch2::Catch2
    fmindex-collection
)


target_include_directories(${PROJECT_NAME} PRIVATE ${DIVSUFSORT64_INCLUDE_DIRS})
target_compile_options(${PROJECT_NAME} PRIVATE ${DIVSUFSORT64_CFLAGS})
target_link_libraries(${PROJECT_NAME} ${DIVSUFSORT64_LIBRARIES})
add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
