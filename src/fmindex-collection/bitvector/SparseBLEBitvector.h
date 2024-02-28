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

namespace fmindex_collection::bitvector {
/**
 * SparseBLEBitvector sparse block length encoded bitvector.
 *
 * Uses Block length encoding and simplification on the
 * assumption the bit vector is sparse
 */
template <int64_t _BlockLengthE = 2, BitVector_c BV1 = Bitvector, BitVector_c BV2 = Bitvector>
struct SparseBLEBitvector {

    static constexpr bool   CompressOnes = _BlockLengthE<0;
    static constexpr size_t BlockLengthE = CompressOnes?-_BlockLengthE:_BlockLengthE;
    BV1 indicatorBitvector;
    BV2 uncompressedBitvector;
    size_t totalLength{};

    std::vector<uint8_t> trailing; // extra characters, that aren't inside a block

    template <typename CB>
    SparseBLEBitvector(size_t length, CB cb)
        : SparseBLEBitvector{std::views::iota(size_t{}, length) | std::views::transform([&](size_t i) {
            return cb(i);
        })}
    {}

    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint8_t>
    SparseBLEBitvector(range_t&& _range) {
        for (auto v : _range) {
            push_back(v);
        }
    }

    SparseBLEBitvector() = default;
    SparseBLEBitvector(SparseBLEBitvector const&) = default;
    SparseBLEBitvector(SparseBLEBitvector&&) noexcept = default;
    auto operator=(SparseBLEBitvector const&) -> SparseBLEBitvector& = default;
    auto operator=(SparseBLEBitvector&&) noexcept -> SparseBLEBitvector& = default;

    void push_back(bool _value) {
        trailing.push_back(_value);
        if (trailing.size() == (1<<BlockLengthE)) {
            bool allTheSame = true;
            bool first = trailing[0];
            for (auto v : trailing) {
                allTheSame &= (first == v);
            }
            if constexpr (!CompressOnes) {
                if (first) allTheSame = false; // Don't compress ones
            } else {
                if (!first) allTheSame = false; // Don't compres zeros
            }
            indicatorBitvector.push_back(allTheSame);
            if (!allTheSame) {
                for (size_t j{0}; j < (1<<BlockLengthE); ++j) {
                    uncompressedBitvector.push_back(trailing[j]);
                }
            }
            totalLength += trailing.size();
            trailing.clear();
        }
    }

    size_t size() const noexcept {
        return totalLength + trailing.size();
    }

    bool symbol(size_t idx) const noexcept {
        if (idx >= totalLength) {
            return trailing[idx-totalLength];
        }
        auto blockId = idx >> BlockLengthE;

        // If compressed, must be what ever is compressed
        if (indicatorBitvector.symbol(blockId)) {
            return CompressOnes;
        }

        auto compBlocks = indicatorBitvector.rank(blockId);

        auto detailSymbId = idx - (compBlocks * (1 << BlockLengthE));
        return uncompressedBitvector.symbol(detailSymbId);
    }

    uint64_t rank(size_t idx) const noexcept {
        if (idx >= totalLength) {
            auto r = uncompressedBitvector.rank(uncompressedBitvector.size());
            if constexpr (CompressOnes) {
                auto compBlocks = indicatorBitvector.rank(indicatorBitvector.size());
                r += (compBlocks * (1 << BlockLengthE));
            }
            for (size_t i{0}; i < idx-totalLength; ++i) {
                r += trailing[i];
            }
            return r;
        }

        auto blockId        = idx >> BlockLengthE;

        // If target block is compressed
        if (indicatorBitvector.symbol(blockId)) {
            auto compBlocks     = indicatorBitvector.rank(blockId);
            auto detailedSymbId = ((idx >> BlockLengthE) - compBlocks) << BlockLengthE;
            auto r = uncompressedBitvector.rank(detailedSymbId);

            if constexpr (CompressOnes) {
                r += (compBlocks * (1 << BlockLengthE));
                r += idx % (1 << BlockLengthE);
            }
            return r;
        } else {
            auto compBlocks     = indicatorBitvector.rank(blockId);
            auto detailedSymbId = idx - (compBlocks * (1 << BlockLengthE));
            auto r = uncompressedBitvector.rank(detailedSymbId);
            if constexpr (CompressOnes) {
                r += (compBlocks * (1 << BlockLengthE));
            }
            return r;
        }
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(indicatorBitvector, uncompressedBitvector, totalLength, trailing);
    }
};
static_assert(BitVector_c<SparseBLEBitvector<>>);

}
