// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "concepts.h"

#include <array>
#include <bitset>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <ranges>
#include <span>
#include <vector>

namespace fmindex_collection::bitvector {

/**
 * L1Bitvector a bit vector with only bits and blocks
 *
 *   For 64bits, we need 128bits, resulting in 2bits per bit
 */
struct L1Bitvector {
    std::vector<uint64_t> superblocks{0, 0};
    std::vector<uint64_t> bits{0};
    size_t totalLength{};

    L1Bitvector() = default;
    L1Bitvector(L1Bitvector const&) = default;
    L1Bitvector(L1Bitvector&&) noexcept = default;

    template <typename CB>
    L1Bitvector(size_t length, CB cb)
        : L1Bitvector{std::views::iota(size_t{}, length) | std::views::transform([&](size_t i) {
            return cb(i);
        })}
    {}

    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint8_t>
    L1Bitvector(range_t&& _range) {
        reserve(_range.size());

        auto iter = _range.begin();

        size_t const loop64  = _range.size() / 64;
        for (size_t l64{}; l64 < loop64; ++l64) {
            for (size_t i{}; i < 64; ++i) {
                bool value         = *(iter++);
                bits.back()        |= (uint64_t{value} << i);
                superblocks.back() += value;
            }
            totalLength += 64;
            superblocks.emplace_back(superblocks.back());
            bits.emplace_back();
        }
        while(totalLength < _range.size()) {
            push_back(*(iter++));
        }
    }

    auto operator=(L1Bitvector const&) -> L1Bitvector& = default;
    auto operator=(L1Bitvector&&) noexcept -> L1Bitvector& = default;

    void reserve(size_t _length) {
        superblocks.reserve((_length+1)/64 + 1);
        bits.reserve(_length/64 + 1);
    }

    void push_back(bool _value) {
        if (_value) {
            auto bitId         = totalLength % 64;
            bits.back()        |= (uint64_t{_value} << bitId);
            superblocks.back() += _value;
        }
        totalLength += 1;
        if (totalLength % 64 == 0) { // new superblock
            superblocks.emplace_back(superblocks.back());
            bits.emplace_back();
        }
    }

    size_t size() const noexcept {
        return totalLength;
    }

    bool symbol(size_t idx) const noexcept {
        auto bitId        = idx % 64;
        auto blockId      = idx / 64;
        auto bit = (bits[blockId] >> bitId) & 1;
        return bit;
    }

    uint64_t rank(size_t idx) const noexcept {
        auto bitId        = idx % 64;
        auto superblockId = idx / 64;
        if (idx % 64 == 0) return superblocks[superblockId];

        auto maskedBits = (bits[superblockId] << (64 - bitId));
        auto bitcount   = std::bitset<64>{maskedBits}.count();

        return superblocks[superblockId]
                + bitcount;
    }



    template <typename Archive>
    void serialize(Archive& ar) {
        ar(superblocks, bits, totalLength);
    }
};

static_assert(BitVector_c<L1Bitvector>);

}
