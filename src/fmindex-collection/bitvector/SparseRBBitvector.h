// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "concepts.h"
#include "Bitvector.h"

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
 * SparseRBBitvector sparse block length encoded bitvector.
 *
 * Uses Block length encoding and simplification on the
 * assumption the bit vector is sparse
 */
template <int64_t BlockLengthE = 2, Bitvector_c BV1 = Bitvector, Bitvector_c BV2 = Bitvector>
struct SparseRBBitvector {

    static constexpr size_t BlockLength = (size_t{1} << BlockLengthE);
    BV1 indicatorBitvector;
    BV2 uncompressedBitvector;
    size_t totalLength{};

    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint8_t>
    SparseRBBitvector(range_t&& _range) {
        size_t completeBlocks = _range.size() / BlockLength;

        // Save data in temporary structure
        auto tmpIndicatorBV = std::vector<bool>{};
        auto tmpUncompressedBV = std::vector<bool>{};

        for (size_t block{0}; block < completeBlocks; ++block) {
            bool allZeros = true;
            for (size_t j{0}; j < BlockLength; ++j) {
                allZeros &= (0 == _range[block*BlockLength + j]);
            }
            tmpIndicatorBV.push_back(allZeros);
            if (!allZeros) {
                for (size_t j{0}; j < BlockLength; ++j) {
                    tmpUncompressedBV.push_back(_range[block*BlockLength+j]);
                }
            }
        }

        // last block is not compressed, need a last block, even if empty.
        // required to allow rank(size()) access
        tmpIndicatorBV.push_back(0);
        for (size_t j{0}; j < BlockLength; ++j) {
            auto k = completeBlocks*BlockLength+j;
            if (k < _range.size()) {
                tmpUncompressedBV.push_back(_range[k]);
            } else {
                tmpUncompressedBV.push_back(0);
            }
        }

        indicatorBitvector = BV1{tmpIndicatorBV};
        uncompressedBitvector = BV2{tmpUncompressedBV};
        totalLength = _range.size();
    }

    SparseRBBitvector() = default;
    SparseRBBitvector(SparseRBBitvector const&) = default;
    SparseRBBitvector(SparseRBBitvector&&) noexcept = default;

    auto operator=(SparseRBBitvector const&) -> SparseRBBitvector& = default;
    auto operator=(SparseRBBitvector&&) noexcept -> SparseRBBitvector& = default;

    size_t size() const noexcept {
        return totalLength;
    }

    bool symbol(size_t idx) const noexcept {
        assert(idx < size());
        auto blockId = idx >> BlockLengthE;

        // If compressed, must be what ever is compressed
        if (indicatorBitvector.symbol(blockId)) {
            return 0;
        }

        auto compBlocks = indicatorBitvector.rank(blockId);

        auto detailSymbId = idx - (compBlocks * BlockLength);
        return uncompressedBitvector.symbol(detailSymbId);
    }

    uint64_t rank(size_t idx) const noexcept {
        assert(idx <= size());

        auto blockId    = idx >> BlockLengthE;
        auto [ind, compBlocks] = [&]() -> std::tuple<bool, uint64_t> {
            auto ind        = indicatorBitvector.symbol(blockId);
            auto compBlocks = indicatorBitvector.rank(blockId);
            return {ind, compBlocks};
        }();

        // Remove bits in last block, if block is compressed
        auto detailedSymbId = idx - (idx % BlockLength) * ind - (compBlocks * BlockLength);
        auto r = uncompressedBitvector.rank(detailedSymbId);
        assert(r <= idx);
        return r;
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(indicatorBitvector, uncompressedBitvector, totalLength);
    }

    static size_t estimateSize(size_t blockCt, size_t zeroBlocks, size_t oneBlocks) {
        (void)oneBlocks; // unused, not required for estimation
        auto nonZeroBlocks = blockCt - zeroBlocks;
        return BV1::estimateSize(blockCt) + BV2::estimateSize(nonZeroBlocks*BlockLength);
    }
};
static_assert(Bitvector_c<SparseRBBitvector<>>);

}
