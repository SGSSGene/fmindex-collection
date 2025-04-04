// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../bitset_popcount.h"
#include "../utils.h"
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
 * PairedL0_NBitvector a bit vector with only bits and blocks
 *
 *   (64) For 128bits, we need 192bits, resulting in 1.5bits per bit
 *   (128) for 256bits, we need 320bits, resulting in 1.25bits per bit
 *   (256) for 512bits, we need 576bits, resulting in 1.125bits per bit
 *   (512) for 1024bits, we need 1088bits, resulting in 1.0625bits per bit
 *   (1024) for 2048bits, we need 2112bits, resulting in 1.0312bits per bit
 */
template <size_t bits_ct, bool Align=true>
struct PairedL0_NBitvector {
    std::vector<uint64_t>                      l0{0};
    std::vector<AlignedBitset<bits_ct, Align>> bits{{}};
    size_t totalLength{};

    PairedL0_NBitvector() = default;
    PairedL0_NBitvector(PairedL0_NBitvector const&) = default;
    PairedL0_NBitvector(PairedL0_NBitvector&&) noexcept = default;

    template <typename CB>
    PairedL0_NBitvector(size_t length, CB cb)
        : PairedL0_NBitvector{std::views::iota(size_t{}, length) | std::views::transform([&](size_t i) {
            return cb(i);
        })}
    {}

    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint8_t>
    PairedL0_NBitvector(range_t&& _range) {
        reserve(_range.size());

        auto iter = _range.begin();

        size_t const loop64  = _range.size() / (bits_ct*2);

        for (size_t l64{}; l64 < loop64; ++l64) {
            for (size_t i{}; i < bits_ct; ++i) {
                bool value         = *(iter++);
                bits.back()[i]     = value;
                l0.back() += value;
            }
            bits.emplace_back();
            l0.emplace_back(l0.back());
            for (size_t i{}; i < bits_ct; ++i) {
                bool value         = *(iter++);
                bits.back()[i]     = value;
                l0.back() += value;
            }
            bits.emplace_back();
        }
        totalLength = bits_ct * (loop64*2);
        while(totalLength < _range.size()) {
            push_back(*(iter++));
        }
    }

    auto operator=(PairedL0_NBitvector const&) -> PairedL0_NBitvector& = default;
    auto operator=(PairedL0_NBitvector&&) noexcept -> PairedL0_NBitvector& = default;

    void reserve(size_t _length) {
        l0.reserve((_length+1)/(bits_ct*2) + 1);
        bits.reserve(_length/bits_ct + 1);
    }

    void push_back(bool _value) {
        if (_value) {
            auto bitId         = totalLength % bits_ct;
            bits.back()[bitId] = _value;
            l0.back() += _value;
        }
        totalLength += 1;
        if (totalLength % (bits_ct*2) == 0) { // new superblock
            bits.emplace_back();
        } else if (totalLength % (bits_ct*2) == bits_ct) { // new superblock
            l0.emplace_back(l0.back());
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
        auto bitId = idx % (bits_ct*2);
        auto superblockId = idx / bits_ct;

        auto right = (superblockId%2);
        auto count = signed_rshift_and_count(bits[superblockId].bits, bitId);

        // Implicit conversions, because emcc can't handle over/underflow correctly
        auto ct = static_cast<int64_t>(l0[superblockId/2]) + (static_cast<int64_t>(right)*2-1) * static_cast<int64_t>(count);
        return ct;
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(l0, totalLength, bits);
    }
};

using PairedL0_64Bitvector  = PairedL0_NBitvector<64>;
using PairedL0_128Bitvector = PairedL0_NBitvector<128>;
using PairedL0_256Bitvector = PairedL0_NBitvector<256>;
using PairedL0_512Bitvector = PairedL0_NBitvector<512>;
using PairedL0_1024Bitvector = PairedL0_NBitvector<1024>;
using PairedL0_2048Bitvector = PairedL0_NBitvector<2048>;

static_assert(BitVector_c<PairedL0_64Bitvector>);
static_assert(BitVector_c<PairedL0_128Bitvector>);
static_assert(BitVector_c<PairedL0_256Bitvector>);
static_assert(BitVector_c<PairedL0_512Bitvector>);
static_assert(BitVector_c<PairedL0_1024Bitvector>);
static_assert(BitVector_c<PairedL0_2048Bitvector>);

using PairedL0_64BitvectorUA  = PairedL0_NBitvector<64, false>;
using PairedL0_128BitvectorUA = PairedL0_NBitvector<128, false>;
using PairedL0_256BitvectorUA = PairedL0_NBitvector<256, false>;
using PairedL0_512BitvectorUA = PairedL0_NBitvector<512, false>;
using PairedL0_1024BitvectorUA = PairedL0_NBitvector<1024, false>;
using PairedL0_2048BitvectorUA = PairedL0_NBitvector<2048, false>;

static_assert(BitVector_c<PairedL0_64BitvectorUA>);
static_assert(BitVector_c<PairedL0_128BitvectorUA>);
static_assert(BitVector_c<PairedL0_256BitvectorUA>);
static_assert(BitVector_c<PairedL0_512BitvectorUA>);
static_assert(BitVector_c<PairedL0_1024BitvectorUA>);
static_assert(BitVector_c<PairedL0_2048BitvectorUA>);

}
