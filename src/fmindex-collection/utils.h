#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <sdsl/divsufsort.hpp>
#include <stdexcept>
#include <vector>

namespace fmindex_collection {

inline auto createSA(std::vector<uint8_t> const& input) -> std::vector<int64_t> {
    auto sa = std::vector<int64_t>{};
    sa.resize(input.size());
    auto error = sdsl::divsufsort<int64_t>(static_cast<uint8_t const*>(input.data()), sa.data(), input.size());
    if (error != 0) {
        throw std::runtime_error("some error while creating the suffix array");
    }
    return sa;
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


}
