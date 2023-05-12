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
#include <libsais64.h>
#include <numeric>
#include <span>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace fmindex_collection {

inline auto createSA(std::span<uint8_t const> input, size_t threadNbr) -> std::vector<int64_t> {
    auto sa = std::vector<int64_t>(input.size());
    if (input.size() == 0) {
        return sa;
    }
#if LIBSAIS_OPENMP
    auto r = libsais64_omp(input.data(), sa.data(), input.size(), 0, nullptr, threadNbr);
#else
    (void)threadNbr; // Unused if no openmp is available
    auto r = libsais64(input.data(), sa.data(), input.size(), 0, nullptr);
#endif

    if (r != 0) { throw std::runtime_error("something went wrong constructing the SA"); }
    return sa;
}


inline auto createBWT(std::span<uint8_t const> input, std::span<int64_t const> sa) -> std::vector<uint8_t> {
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
        auto textLen    = l.size();
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
        // number of delimiters ('$') which need to be added. It must be at least one, and it
        // has to make sure the text will be a multiple of samplingRate
        size_t delimCount = samplingRate - (ls % samplingRate);
        inputText.resize(inputText.size() + ls + delimCount, 0);

        if (not reverse) {
            std::ranges::copy(l, end(inputText) - ls - delimCount);
        } else {
//!TODO hack for clang, broken in clang 15
#if __clang__
            auto l2 = std::vector<uint8_t>(l);
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
