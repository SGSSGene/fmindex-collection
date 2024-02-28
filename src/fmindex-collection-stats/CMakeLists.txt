# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0
cmake_minimum_required (VERSION 3.25)


project(fmindex-collection-stats LANGUAGES CXX
        DESCRIPTION "Generates statistics about different data structures.")

add_executable(${PROJECT_NAME}
    main.cpp
)
target_link_libraries(${PROJECT_NAME}
    PRIVATE
    fmindex-collection::fmindex-collection
    nanobench::nanobench
    fmt::fmt-header-only
)
