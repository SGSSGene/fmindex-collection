# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0
cmake_minimum_required (VERSION 3.8)

project(test_search_schemes)

# search scheme generator unittests
add_executable(${PROJECT_NAME}
    checkGenerators.cpp
    checkGeneratorsIsComplete.cpp
    expand.cpp
    isComplete.cpp
    isValid.cpp
    nodeCount.cpp
    weightedNodeCount.cpp
)
target_link_libraries(${PROJECT_NAME}
    Catch2::Catch2WithMain
    search_schemes
)
add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
