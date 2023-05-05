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

auto createSequences(Sequences auto const& _input, int samplingRate, bool reverse=false) -> std::tuple<size_t, std::vector<uint8_t>, std::vector<size_t>> {
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
    auto inputSizes = std::vector<size_t>{};
    inputSizes.reserve(_input.size());


    for (auto const& l : _input) {
        if (not reverse) {
            inputText.resize(inputText.size() + l.size());
            std::copy_backward(begin(l), end(l), end(inputText));
        } else {
//!TODO hack for clang, broken in clang 15
#if __clang__
            auto l2 = l;
            std::ranges::reverse(l2);
#else
            auto l2 = std::views::reverse(l);
#endif
            inputText.insert(end(inputText), begin(l2), end(l2));
        }
        inputSizes.emplace_back(l.size()+1);
        inputText.emplace_back(0);

        while (inputText.size() % samplingRate) {
            inputSizes.back() += 1;
            inputText.emplace_back(0);
        }
    }
    return {totalSize, inputText, inputSizes};
}

}
