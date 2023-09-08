// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "concepts.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <libsais.h>
#include <libsais64.h>
#include <numeric>
#include <span>
#include <stdexcept>
#include <tuple>
#include <vector>

#if LIBSAIS_OPENMP
#   include <omp.h>
#endif

namespace fmindex_collection {

inline auto createSA(std::span<uint8_t const> input, size_t threadNbr) -> std::vector<uint64_t> {
    auto sa = std::vector<uint64_t>(input.size());
    if (input.size() == 0) {
        return sa;
    }
#if LIBSAIS_OPENMP
    auto r = libsais64_omp(input.data(), reinterpret_cast<int64_t*>(sa.data()), input.size(), 0, nullptr, threadNbr);
#else
    (void)threadNbr; // Unused if no openmp is available
    auto r = libsais64(input.data(), reinterpret_cast<int64_t*>(sa.data()), input.size(), 0, nullptr);
#endif

    if (r != 0) { throw std::runtime_error("something went wrong constructing the SA"); }
    return sa;
}


inline auto createBWT(std::span<uint8_t const> input, std::span<uint64_t const> sa) -> std::vector<uint8_t> {
    assert(input.size() == sa.size());
    auto bwt = std::vector<uint8_t>{};
    bwt.resize(input.size());
    for (size_t i{0}; i < sa.size(); ++i) {
        bwt[i] = input[(sa[i] + input.size() - 1) % input.size()];
    }
    return bwt;
}

auto createSequences(Sequences auto const& _input, int samplingRate, bool reverse=false) -> std::tuple<size_t, std::vector<uint8_t>, std::vector<std::tuple<size_t, size_t>>> {
    // compute total numbers of bytes of the text including delimiters "$"
    size_t totalSize{};
    for (auto const& l : _input) {
        auto textLen  = l.size();
        auto delimLen = samplingRate - textLen % samplingRate; // Make sure it is always a multiple of samplingRate
        totalSize += textLen + delimLen;
    }

    // our concatenated sequences with delimiters
    auto inputText = std::vector<uint8_t>{};
    inputText.reserve(totalSize);

    // list of sizes of the individual sequences
    auto inputSizes = std::vector<std::tuple<size_t, size_t>>{};
    inputSizes.reserve(_input.size());

    for (auto const& l : _input) {
        auto ls = l.size();

        if (not reverse) {
            inputText.insert(inputText.end(), begin(l), end(l));
        } else {
//!TODO hack for clang, broken in clang 15
#if __clang__
            auto l2 = std::vector<uint8_t>(l);
            std::ranges::reverse(l2);
#else
            auto l2 = std::views::reverse(l);
#endif
            inputText.insert(inputText.end(), begin(l2), end(l2));
        }

        // number of delimiters ('$') which need to be added. It must be at least one, and it
        // has to make sure the text will be a multiple of samplingRate
        size_t delimCount = samplingRate - (ls % samplingRate);

        // fill with delimiters/zeros
        inputText.resize(inputText.size() + delimCount);

        inputSizes.emplace_back(ls, delimCount);
    }
    return {totalSize, inputText, inputSizes};
}

auto createSequencesAndReverse(Sequences auto const& _input, int samplingRate) -> std::tuple<size_t, std::vector<uint8_t>, std::vector<std::tuple<size_t, size_t>>> {
    // compute total numbers of bytes of the text including delimiters "$"
    size_t totalSize{};
    for (auto const& l : _input) {
        auto textLen  = l.size();
        auto delimLen = samplingRate - textLen % samplingRate; // Make sure it is always a multiple of samplingRate
        totalSize += textLen + delimLen;
    }
    for (auto const& l : std::views::reverse(_input)) {
        auto textLen  = l.size();
        auto delimLen = samplingRate - textLen % samplingRate; // Make sure it is always a multiple of samplingRate
        totalSize += textLen + delimLen;
    }

    // our concatenated sequences with delimiters
    auto inputText = std::vector<uint8_t>{};
    inputText.reserve(totalSize);

    // list of sizes of the individual sequences
    auto inputSizes = std::vector<std::tuple<size_t, size_t>>{};
    inputSizes.reserve(_input.size());

    // add text
    for (auto const& l : _input) {
        auto ls = l.size();
        inputText.insert(inputText.end(), begin(l), end(l));

        // number of delimiters ('$') which need to be added. It must be at least one, and it
        // has to make sure the text will be a multiple of samplingRate
        size_t delimCount = samplingRate - (ls % samplingRate);

        // fill with delimiters/zeros
        inputText.resize(inputText.size() + delimCount);

        inputSizes.emplace_back(ls, delimCount);
    }

    // add reversed text
    for (auto const& l : std::views::reverse(_input)) {
        auto ls = l.size();
        auto l2 = std::views::reverse(l);
        inputText.insert(inputText.end(), begin(l2), end(l2));

        // number of delimiters ('$') which need to be added. It must be at least one, and it
        // has to make sure the text will be a multiple of samplingRate
        size_t delimCount = samplingRate - (ls % samplingRate);

        // fill with delimiters/zeros
        inputText.resize(inputText.size() + delimCount);

        inputSizes.emplace_back(ls, delimCount);
    }

    return {totalSize, inputText, inputSizes};
}



inline auto createSA_32(std::span<uint32_t const> input, size_t threadNbr) -> std::vector<int32_t> {
    auto sa = std::vector<int32_t>(input.size());
    if (input.size() == 0) {
        return sa;
    }
#if LIBSAIS_OPENMP
    auto r = libsais_int_omp((int32_t*)input.data(), sa.data(), input.size(), 0, nullptr, threadNbr);
#else
    (void)threadNbr; // Unused if no openmp is available
    auto r = libsais_int((int32_t*)input.data(), sa.data(), input.size(), 65536, 0);
#endif

    if (r != 0) { throw std::runtime_error("something went wrong constructing the SA"); }
    return sa;
}


inline auto createBWT_32(std::span<uint32_t const> input, std::span<int32_t const> sa) -> std::vector<uint32_t> {
    assert(input.size() == sa.size());
    auto bwt = std::vector<uint32_t>{};
    bwt.resize(input.size());
    for (size_t i{0}; i < sa.size(); ++i) {
        bwt[i] = input[(sa[i] + input.size() - 1) % input.size()];
    }
    return bwt;
}

auto createSequences_32(Sequences auto const& _input, int samplingRate, bool reverse=false) -> std::tuple<size_t, std::vector<uint32_t>, std::vector<std::tuple<size_t, size_t>>> {
    // compute total numbers of bytes of the text including delimiters "$"
    size_t totalSize{};
    for (auto const& l : _input) {
        auto textLen    = l.size();
        auto delimLen = samplingRate - textLen % samplingRate; // Make sure it is always a multiple of samplingRate
        totalSize += textLen + delimLen;
    }

    // our concatenated sequences with delimiters
    auto inputText = std::vector<uint32_t>{};
    inputText.reserve(totalSize);

    // list of sizes of the individual sequences
    auto inputSizes = std::vector<std::tuple<size_t, size_t>>{};
    inputSizes.reserve(_input.size());


    for (auto const& l : _input) {
        auto ls = l.size();
        // number of delimiters ('$') which need to be added. It must be at least one, and it
        // has to make sure the text will be a multiple of samplingRate
        size_t delimCount = samplingRate - (ls % samplingRate);
        inputText.resize(inputText.size() + ls + delimCount, 0);

        if (not reverse) {
            std::ranges::copy(l, end(inputText) - ls - delimCount);
        } else {
//!TODO hack for clang, broken in clang 15
#if __clang__
            auto l2 = std::vector<uint32_t>(l);
            std::ranges::reverse(l2);
#else
            auto l2 = std::views::reverse(l);
#endif
            std::ranges::copy(l2, end(inputText) - ls - delimCount);
        }

        inputSizes.emplace_back(l.size(), delimCount);
    }
    return {totalSize, inputText, inputSizes};
}


}
