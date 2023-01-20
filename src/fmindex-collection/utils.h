#pragma once

#include "concepts.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <libsaispp/sais.hpp>
#include <numeric>
#include <span>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace fmindex_collection {

inline auto createSA(std::vector<uint8_t> const& input) -> std::vector<int64_t> {
    return libsais::constructSA64(input);
}


inline auto createBWT(std::vector<uint8_t> const& input, std::vector<int64_t> const& sa) -> std::vector<uint8_t> {
    assert(input.size() == sa.size());
    std::vector<uint8_t> bwt;
    bwt.resize(input.size());
    for (size_t i{0}; i < sa.size(); ++i) {
        bwt[i] = input[(sa[i] + input.size()- 1) % input.size()];
    }
    return bwt;
}

auto createSequences(Sequences auto const& _input, bool reverse=false) -> std::tuple<size_t, std::vector<uint8_t>, std::vector<size_t>> {
    // compute total numbers of bytes of the text including delimiters "$"
    size_t totalSize = std::accumulate(begin(_input), end(_input), size_t{0}, [](auto s, auto const& l) { return s + l.size() + 1; });

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
        inputText.emplace_back(0);
        inputSizes.emplace_back(l.size()+1);
    }
    return {totalSize, inputText, inputSizes};
}

}
