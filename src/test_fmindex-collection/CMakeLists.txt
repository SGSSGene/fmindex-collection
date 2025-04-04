# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0
cmake_minimum_required (VERSION 3.10)

project(test_fmindex-collection)

add_executable(${PROJECT_NAME}
    bitvector/bitset_popcount.cpp
    bitvector/benchmark.cpp
    bitvector/unittest.cpp
    checkTernarylogic.cpp
    fmindex/checkBiFMIndex.cpp
    fmindex/checkBiFMIndexCursor.cpp
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
    search/checkLocateFMTree.cpp
    search/checkReverseIndexSearch.cpp
    search/checkSearchBacktracking.cpp
    search/checkSearchPseudo.cpp
    search/checkSearches.cpp
    search/checkSearchHammingSM.cpp
    search_scheme/checkGenerators.cpp
    search_scheme/checkGeneratorsIsComplete.cpp
    search_scheme/expand.cpp
    search_scheme/isComplete.cpp
    search_scheme/isValid.cpp
    search_scheme/nodeCount.cpp
    search_scheme/weightedNodeCount.cpp
    string/benchmark_alphabet_4.cpp
    string/benchmark_alphabet_5.cpp
    string/benchmark_alphabet_6.cpp
    string/benchmark_alphabet_8.cpp
    string/benchmark_alphabet_16.cpp
    string/benchmark_alphabet_32.cpp
    string/benchmark_alphabet_64.cpp
    string/benchmark_alphabet_128.cpp
    string/benchmark_alphabet_255.cpp
    string/unittest.cpp
    suffixarray/checkCSA.cpp
    utils.cpp
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
if(FMC_USE_SUX AND FMC_USE_PASTA AND FMC_USE_SDSL)
    # bitvector/benchmark requires sux and paste to compile successful
    target_compile_definitions(${PROJECT_NAME} INTERFACE -DFMC_USE_SUX)
    target_compile_definitions(${PROJECT_NAME} INTERFACE -DFMC_USE_PASTA)
    target_sources(${PROJECT_NAME} PRIVATE bitvector/benchmark.cpp)
    target_link_libraries(${PROJECT_NAME}
        pasta_bit_vector
        sux::sux
    )
endif()

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
