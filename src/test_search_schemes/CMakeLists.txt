# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file.
# -----------------------------------------------------------------------------------------------------
cmake_minimum_required (VERSION 3.8)

project(test_search_schemes)

# search scheme generator unittests
add_executable(${PROJECT_NAME}
    checkGenerators.cpp
    checkGeneratorsIsComplete.cpp
    expand.cpp
    isComplete.cpp
    isValid.cpp
    main.cpp
    nodeCount.cpp
)
target_link_libraries(${PROJECT_NAME}
    Catch2::Catch2
    search_schemes
)
add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
