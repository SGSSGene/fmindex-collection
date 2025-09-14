// SPDX-FileCopyrightText: 2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2025 Knut Reinert & MPI für molekulare Genetik
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
 * RBBitvector sparse block length encoded bitvector.
 *
 * Similar to SparseBLEBitvector, but compresses blocks of 0 and 1
 */
template <int64_t BlockLengthE = 2, Bitvector_c BV1 = Bitvector, Bitvector_c BV2 = Bitvector>
struct RBBitvector {

    static constexpr size_t BlockLength = (size_t{1} << BlockLengthE);
    BV1 indicatorBitvector;
    BV2 uncompressedBitvector;
    BV2 zerosOrOnesBitvector;
    size_t totalLength{};

    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint8_t>
    RBBitvector(range_t&& _range) {
        size_t completeBlocks = _range.size() / BlockLength;

        // Save data in temporary structure
        auto tmpIndicatorBV = std::vector<bool>{};
        auto tmpUncompressedBV = std::vector<bool>{};
        auto tmpZerosOrOnesBV = std::vector<bool>{};

        for (size_t block{0}; block < completeBlocks; ++block) {
            bool allZeros = true;
            bool allOnes = true;
            for (size_t j{0}; j < BlockLength; ++j) {
                allZeros &= (0 == _range[block*BlockLength + j]);
                allOnes &= (1 == _range[block*BlockLength + j]);
            }
            bool ind = allZeros || allOnes;
            tmpIndicatorBV.push_back(ind);
            if (!ind) {
                for (size_t j{0}; j < BlockLength; ++j) {
                    tmpUncompressedBV.push_back(_range[block*BlockLength+j]);
                }
            } else {
                tmpZerosOrOnesBV.push_back(allOnes);
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
        zerosOrOnesBitvector = BV2{tmpZerosOrOnesBV};
        totalLength = _range.size();
    }

    RBBitvector() = default;
    RBBitvector(RBBitvector const&) = default;
    RBBitvector(RBBitvector&&) noexcept = default;

    auto operator=(RBBitvector const&) -> RBBitvector& = default;
    auto operator=(RBBitvector&&) noexcept -> RBBitvector& = default;

    size_t size() const noexcept {
        return totalLength;
    }

    bool symbol(size_t idx) const noexcept {
        assert(idx < size());
        auto blockId = idx >> BlockLengthE;

        auto ind        = indicatorBitvector.symbol(blockId);
        auto compBlocks = indicatorBitvector.rank(blockId);

        // If compressed, must be what ever is compressed
        if (ind) {
            return zerosOrOnesBitvector.symbol(compBlocks);
        }

        auto detailSymbId = idx - (compBlocks * BlockLength);
        return uncompressedBitvector.symbol(detailSymbId);
    }

    uint64_t rank(size_t idx) const noexcept {
        assert(idx <= size());

        auto blockId    = idx >> BlockLengthE;
        auto ind        = indicatorBitvector.symbol(blockId);
        auto compBlocks = indicatorBitvector.rank(blockId);

        auto compOnes   = zerosOrOnesBitvector.rank(compBlocks) * BlockLength;
        auto detailedSymbId = idx - (compBlocks * BlockLength);
        if (ind == 0) {
            auto uncompOnes = uncompressedBitvector.rank(detailedSymbId);
            auto r = compOnes + uncompOnes;
            assert(r <= idx);
            return r;
        }
        auto uncompOnes = uncompressedBitvector.rank(detailedSymbId - (idx % BlockLength));
        auto allOnes = zerosOrOnesBitvector.symbol(compBlocks);
        auto r = compOnes + uncompOnes + (idx % BlockLength) * allOnes;
        assert(r <= idx);
        return r;
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(indicatorBitvector, uncompressedBitvector, zerosOrOnesBitvector, totalLength);
    }

    static size_t estimateSize(size_t totalSize, size_t zeroBlocks, size_t oneBlocks) {
        size_t blockCt = (totalSize+BlockLength-1)/BlockLength;
        auto mixedBlocks = blockCt - zeroBlocks - oneBlocks;
        return BV1::estimateSize(blockCt) + BV2::estimateSize(zeroBlocks + oneBlocks) + BV2::estimateSize(mixedBlocks*BlockLength);
    }

};
static_assert(Bitvector_c<RBBitvector<>>);

}
