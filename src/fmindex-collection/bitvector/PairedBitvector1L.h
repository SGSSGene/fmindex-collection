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

namespace fmc::bitvector {

/**
 * PairedBitvector1L a bit vector with only bits and blocks
 *
 *   (64) For 128bits, we need 192bits, resulting in 1.5bits per bit
 *   (128) for 256bits, we need 320bits, resulting in 1.25bits per bit
 *   (256) for 512bits, we need 576bits, resulting in 1.125bits per bit
 *   (512) for 1024bits, we need 1088bits, resulting in 1.0625bits per bit
 *   (1024) for 2048bits, we need 2112bits, resulting in 1.0312bits per bit
 */
template <size_t bits_ct, bool Align=true>
struct PairedBitvector1L {
    std::vector<uint64_t>                      l0{0};
    std::vector<AlignedBitset<bits_ct, Align>> bits{{}};
    size_t totalLength{};

    PairedBitvector1L() = default;
    PairedBitvector1L(PairedBitvector1L const&) = default;
    PairedBitvector1L(PairedBitvector1L&&) noexcept = default;

    // constructor accepting view to bools or already compact uint64_t
    template <std::ranges::sized_range range_t>
    PairedBitvector1L(range_t&& _range) {
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
    PairedBitvector1L(range_t&& _range) {
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

    auto operator=(PairedBitvector1L const&) -> PairedBitvector1L& = default;
    auto operator=(PairedBitvector1L&&) noexcept -> PairedBitvector1L& = default;

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
    void serialize(this auto&& self, Archive& ar) {
        ar(self.l0, self.totalLength, self.bits);
    }
};

using PairedBitvector1L_64  = PairedBitvector1L<64>;
using PairedBitvector1L_128 = PairedBitvector1L<128>;
using PairedBitvector1L_256 = PairedBitvector1L<256>;
using PairedBitvector1L_512 = PairedBitvector1L<512>;
using PairedBitvector1L_1024 = PairedBitvector1L<1024>;
using PairedBitvector1L_2048 = PairedBitvector1L<2048>;

static_assert(Bitvector_c<PairedBitvector1L_64>);
static_assert(Bitvector_c<PairedBitvector1L_128>);
static_assert(Bitvector_c<PairedBitvector1L_256>);
static_assert(Bitvector_c<PairedBitvector1L_512>);
static_assert(Bitvector_c<PairedBitvector1L_1024>);
static_assert(Bitvector_c<PairedBitvector1L_2048>);

using PairedBitvector1L_64UA  = PairedBitvector1L<64, false>;
using PairedBitvector1L_128UA = PairedBitvector1L<128, false>;
using PairedBitvector1L_256UA = PairedBitvector1L<256, false>;
using PairedBitvector1L_512UA = PairedBitvector1L<512, false>;
using PairedBitvector1L_1024UA = PairedBitvector1L<1024, false>;
using PairedBitvector1L_2048UA = PairedBitvector1L<2048, false>;

static_assert(Bitvector_c<PairedBitvector1L_64UA>);
static_assert(Bitvector_c<PairedBitvector1L_128UA>);
static_assert(Bitvector_c<PairedBitvector1L_256UA>);
static_assert(Bitvector_c<PairedBitvector1L_512UA>);
static_assert(Bitvector_c<PairedBitvector1L_1024UA>);
static_assert(Bitvector_c<PairedBitvector1L_2048UA>);

}
