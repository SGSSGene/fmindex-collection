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

    // constructor accepting view to bools or already compact uint64_t
    template <std::ranges::sized_range range_t>
    PairedL0_NBitvector(range_t&& _range) {
        auto _size = _range.size();
        if constexpr (std::same_as<std::ranges::range_value_t<range_t>, uint64_t>) {
            *this = {std::forward<range_t>(_range) | view_as_bitset<bits_ct>};
            totalLength = _size*64;
        } else if constexpr (std::convertible_to<std::ranges::range_value_t<range_t>, bool>) {
            *this = {std::forward<range_t>(_range) | view_bool_as_uint64 | view_as_bitset<bits_ct>};
            totalLength = _size;
        } else {
            []<bool b=false>() {
                static_assert(b, "Must be an uint64_t or convertible to bool");
            }();
        }
        l0.resize((totalLength+bits_ct)/(bits_ct*2) + 1);
        bits.resize(totalLength/bits_ct + 1);
    }

    // the actual constructor, already receiving premade std::bitsets<N>
    template <std::ranges::sized_range range_t>
        requires std::same_as<std::ranges::range_value_t<range_t>, std::bitset<bits_ct>>
    PairedL0_NBitvector(range_t&& _range) {
        auto _length = _range.size()*bits_ct;

        l0.resize((_length+bits_ct)/(bits_ct*2) + 1);
        bits.resize(_length/bits_ct + 1);


        for (auto const& b : _range) {
            auto bits_id = totalLength / bits_ct;
            auto l0_id = totalLength / (bits_ct*2);
            bits[bits_id].bits = b;
            totalLength += bits_ct;

            if (totalLength % (bits_ct*2) == bits_ct) { // accumulator for next block
                l0[l0_id] += bits[bits_id].count();
                l0[l0_id+1] = l0[l0_id];
            } else {
                l0[l0_id+1] += bits[bits_id].count();
            }
        }
    }

    auto operator=(PairedL0_NBitvector const&) -> PairedL0_NBitvector& = default;
    auto operator=(PairedL0_NBitvector&&) noexcept -> PairedL0_NBitvector& = default;

    void reserve(size_t _length) {
        l0.reserve((_length+bits_ct)/(bits_ct*2) + 1);
        bits.reserve(_length/bits_ct + 1);
    }

    void push_back(bool _value) {
        auto bitId         = totalLength % bits_ct;
        bits.back()[bitId] = _value;
        l0.back()         += _value;

        totalLength += 1;
        if (totalLength % bits_ct == 0) { // filled a bits block
            if (totalLength % (bits_ct*2) == bits_ct) { // accumulator for next block
                l0.emplace_back(l0.back());
            }
            bits.emplace_back();
        }
    }

    size_t size() const noexcept {
        return totalLength;
    }

    bool symbol(size_t idx) const noexcept {
        assert(idx < totalLength);
        auto bitId        = idx % bits_ct;
        auto blockId      = idx / bits_ct;
        auto bit = bits[blockId][bitId];
        return bit;
    }

    uint64_t rank(size_t idx) const noexcept {
        assert(idx <= totalLength);
        auto bitId = idx % (bits_ct*2);
        auto l0Id  = idx / bits_ct;

        int64_t right_l0 = (l0Id%2)*2-1;

        int64_t count = skip_first_or_last_n_bits_and_count(bits[l0Id].bits, bitId);

        auto ct = l0[l0Id/2] + right_l0 * count;
        assert(ct <= idx);
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

static_assert(Bitvector_c<PairedL0_64Bitvector>);
static_assert(Bitvector_c<PairedL0_128Bitvector>);
static_assert(Bitvector_c<PairedL0_256Bitvector>);
static_assert(Bitvector_c<PairedL0_512Bitvector>);
static_assert(Bitvector_c<PairedL0_1024Bitvector>);
static_assert(Bitvector_c<PairedL0_2048Bitvector>);

using PairedL0_64BitvectorUA  = PairedL0_NBitvector<64, false>;
using PairedL0_128BitvectorUA = PairedL0_NBitvector<128, false>;
using PairedL0_256BitvectorUA = PairedL0_NBitvector<256, false>;
using PairedL0_512BitvectorUA = PairedL0_NBitvector<512, false>;
using PairedL0_1024BitvectorUA = PairedL0_NBitvector<1024, false>;
using PairedL0_2048BitvectorUA = PairedL0_NBitvector<2048, false>;

static_assert(Bitvector_c<PairedL0_64BitvectorUA>);
static_assert(Bitvector_c<PairedL0_128BitvectorUA>);
static_assert(Bitvector_c<PairedL0_256BitvectorUA>);
static_assert(Bitvector_c<PairedL0_512BitvectorUA>);
static_assert(Bitvector_c<PairedL0_1024BitvectorUA>);
static_assert(Bitvector_c<PairedL0_2048BitvectorUA>);

}
