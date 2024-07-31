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
 * L1_NBitvector a bit vector with only bits and blocks
 *
 *   For 64bits,  we need 256bits, resulting in 2.0bits per bit
 *   For 128bits, we need 192bits, resulting in 1.5bits per bit
 *   For 256bits, we need 320bits, resulting in 1.25bits per bit
 *
 */
template <size_t bits_ct>
struct L1_NBitvector {
    std::vector<uint64_t> superblocks{0, 0};
    std::vector<std::bitset<bits_ct>> bits{0};
    size_t totalLength{};

    L1_NBitvector() = default;
    L1_NBitvector(L1_NBitvector const&) = default;
    L1_NBitvector(L1_NBitvector&&) noexcept = default;

    template <typename CB>
    L1_NBitvector(size_t length, CB cb)
        : L1_NBitvector{std::views::iota(size_t{}, length) | std::views::transform([&](size_t i) {
            return cb(i);
        })}
    {}

    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint8_t>
    L1_NBitvector(range_t&& _range) {
        reserve(_range.size());

        auto iter = _range.begin();

        size_t const loop64  = _range.size() / bits_ct;

        for (size_t l64{}; l64 < loop64; ++l64) {
            for (size_t i{}; i < bits_ct; ++i) {
                bool value         = *(iter++);
                bits.back()[i]     = value;
                superblocks.back() += value;
            }
            superblocks.emplace_back(superblocks.back());
            bits.emplace_back();
        }
        totalLength = bits_ct * loop64;
        while(totalLength < _range.size()) {
            push_back(*(iter++));
        }
    }

    auto operator=(L1_NBitvector const&) -> L1_NBitvector& = default;
    auto operator=(L1_NBitvector&&) noexcept -> L1_NBitvector& = default;

    void reserve(size_t _length) {
        superblocks.reserve((_length+1)/bits_ct + 1);
        bits.reserve(_length/bits_ct + 1);
    }

    void push_back(bool _value) {
        if (_value) {
            auto bitId         = totalLength % bits_ct;
            bits.back()[bitId] = _value;
            superblocks.back() += _value;
        }
        totalLength += 1;
        if (totalLength % bits_ct == 0) { // new superblock
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

        auto maskedBits = (bits[superblockId] << (bits_ct - bitId));
        auto bitcount   = maskedBits.count();

        return superblocks[superblockId]
                + bitcount;
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(superblocks, bits, totalLength);
    }
};
using L1_64Bitvector  = L1_NBitvector<64>;
using L1_128Bitvector = L1_NBitvector<128>;
using L1_256Bitvector = L1_NBitvector<256>;

static_assert(BitVector_c<L1_64Bitvector>);
static_assert(BitVector_c<L1_128Bitvector>);
static_assert(BitVector_c<L1_256Bitvector>);

}
