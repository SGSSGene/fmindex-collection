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

#if __has_include(<cereal/types/bitset.hpp>)
#include <cereal/types/bitset.hpp>
#endif

namespace fmindex_collection::bitvector {

/**
 * DoubleL1_NBitvector a bit vector with only bits and blocks
 *
 *   (64) For 128bits, we need 192bits, resulting in 1.5bits per bit
 *   (128) for 256bits, we need 320bits, resulting in 1.25bits per bit
 *   (256) for 512bits, we need 576bits, resulting in 1.125bits per bit
 */
template <size_t bits_ct>
struct DoubleL1_NBitvector {
    std::vector<uint64_t> superblocks{0};
    std::vector<std::bitset<bits_ct>> bits{0};
    size_t totalLength{};

    DoubleL1_NBitvector() = default;
    DoubleL1_NBitvector(DoubleL1_NBitvector const&) = default;
    DoubleL1_NBitvector(DoubleL1_NBitvector&&) noexcept = default;

    template <typename CB>
    DoubleL1_NBitvector(size_t length, CB cb)
        : DoubleL1_NBitvector{std::views::iota(size_t{}, length) | std::views::transform([&](size_t i) {
            return cb(i);
        })}
    {}

    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint8_t>
    DoubleL1_NBitvector(range_t&& _range) {
        reserve(_range.size());

        auto iter = _range.begin();

        size_t const loop64  = _range.size() / (bits_ct*2);

        for (size_t l64{}; l64 < loop64; ++l64) {
            for (size_t i{}; i < bits_ct; ++i) {
                bool value         = *(iter++);
                bits.back()[i]     = value;
                superblocks.back() += value;
            }
            bits.emplace_back();
            superblocks.emplace_back(superblocks.back());
            for (size_t i{}; i < bits_ct; ++i) {
                bool value         = *(iter++);
                bits.back()[i]     = value;
                superblocks.back() += value;
            }
            bits.emplace_back();
        }
        totalLength = bits_ct * (loop64*2);
        while(totalLength < _range.size()) {
            push_back(*(iter++));
        }
    }

    auto operator=(DoubleL1_NBitvector const&) -> DoubleL1_NBitvector& = default;
    auto operator=(DoubleL1_NBitvector&&) noexcept -> DoubleL1_NBitvector& = default;

    void reserve(size_t _length) {
        superblocks.reserve((_length+1)/(bits_ct*2) + 1);
        bits.reserve(_length/bits_ct + 1);
    }

    void push_back(bool _value) {
        if (_value) {
            auto bitId         = totalLength % bits_ct;
            bits.back()[bitId] = _value;
            superblocks.back() += _value;
        }
        totalLength += 1;
        if (totalLength % (bits_ct*2) == 0) { // new superblock
            bits.emplace_back();
        } else if (totalLength % (bits_ct*2) == bits_ct) { // new superblock
            superblocks.emplace_back(superblocks.back());
            bits.emplace_back();
        }

    }

    size_t size() const noexcept {
        return totalLength;
    }

    bool symbol(size_t idx) const noexcept {
        auto bitId        = idx % bits_ct;
        auto blockId      = idx / bits_ct;
        auto bit = bits[blockId][bitId];
        return bit;
    }

    uint64_t rank(size_t idx) const noexcept {
        auto bitId        = idx % bits_ct;
        auto superblockId = idx / bits_ct;

        if (superblockId % 2 == 0) {
            auto maskedBits = (bits[superblockId] >> bitId);
            auto bitcount   = maskedBits.count();

            return superblocks[superblockId/2]
                    - bitcount;
        } else {
            auto maskedBits = (bits[superblockId] << (bits_ct - bitId));
            auto bitcount   = maskedBits.count();

            return superblocks[superblockId/2]
                    + bitcount;
        }
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(superblocks, bits, totalLength);
    }
};
using DoubleL1_64Bitvector  = DoubleL1_NBitvector<64>;
using DoubleL1_128Bitvector = DoubleL1_NBitvector<128>;
using DoubleL1_256Bitvector = DoubleL1_NBitvector<256>;

static_assert(BitVector_c<DoubleL1_64Bitvector>);
static_assert(BitVector_c<DoubleL1_128Bitvector>);
static_assert(BitVector_c<DoubleL1_256Bitvector>);

}
