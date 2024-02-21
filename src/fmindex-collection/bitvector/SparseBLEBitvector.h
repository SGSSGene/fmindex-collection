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
template <size_t BlockLengthE = 2, BitVector_c BV1 = Bitvector, BitVector_c BV2 = Bitvector>
struct SparseBLEBitvector {

    BV1 bv1;
    BV2 bv2;
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
            bool allZeros = true;
            for (size_t j{0}; j < (1<<BlockLengthE); ++j) {
                if (trailing[j]) {
                    allZeros = false;
                    break;
                }
            }
            if (allZeros) {
                bv1.push_back(true);
            } else {
                bv1.push_back(false);
                for (size_t j{0}; j < (1<<BlockLengthE); ++j) {
                    bv2.push_back(trailing[j]);
                }
            }
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

        // If compressed, must be a zero
        if (bv1.symbol(blockId)) {
            return 0;
        }

        auto compBlocks = bv1.rank(blockId);

        auto detailSymbId = idx - (compBlocks * (1 << BlockLengthE));
        return bv2.symbol(detailSymbId);
    }

    uint64_t rank(size_t idx) const noexcept {
        if (idx >= totalLength) {
            auto r = bv2.rank(bv2.size());
            for (auto v : trailing) {
                r += v;
            }
            return r;
        }

        auto blockId = idx >> BlockLengthE;

        auto compBlocks     = bv1.rank(blockId);
        auto detailedSymbId = idx - (compBlocks * (1 << BlockLengthE));

        return bv2.rank(detailedSymbId);
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bv1, bv2, totalLength);
    }
};
static_assert(BitVector_c<SparseBLEBitvector<>>);

}
