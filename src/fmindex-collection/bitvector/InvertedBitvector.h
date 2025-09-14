// SPDX-FileCopyrightText: 2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "concepts.h"
#include "SparseRBBitvector.h"

#include <array>
#include <bitset>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ranges>
#include <stdexcept>
#include <span>
#include <vector>

namespace fmc::bitvector {
/**
 * InvertedBitvector stores the inverted values inside a different bitvector.
 * This is usefull if the other bitvector is a compressed bitvector but only optimized for on-bits or off-bits
 */
template <Bitvector_c BV = SparseRBBitvector<>>
struct InvertedBitvector {
    BV bitvector;

    static constexpr size_t BlockLength = BV::BlockLength;

    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint8_t>
    InvertedBitvector(range_t&& _range)
        : bitvector{_range | std::views::transform([](bool b) { return !b; })} {
    }

    InvertedBitvector() = default;
    InvertedBitvector(InvertedBitvector const&) = default;
    InvertedBitvector(InvertedBitvector&&) noexcept = default;

    auto operator=(InvertedBitvector const&) -> InvertedBitvector& = default;
    auto operator=(InvertedBitvector&&) noexcept -> InvertedBitvector& = default;

    size_t size() const noexcept {
        return bitvector.size();
    }

    bool symbol(size_t idx) const noexcept {
        return !bitvector.symbol(idx);
    }

    uint64_t rank(size_t idx) const noexcept {
        return idx - bitvector.rank(idx);
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bitvector);
    }

    static size_t estimateSize(size_t blockCt, size_t zeroBlocks, size_t oneBlocks) {
        return BV::estimateSize(blockCt, oneBlocks, zeroBlocks);
    }
};
static_assert(Bitvector_c<InvertedBitvector<>>);

}
