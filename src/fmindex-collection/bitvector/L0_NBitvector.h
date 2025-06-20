// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

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
 * L0_NBitvector a bit vector with only bits and blocks
 *
 *   For 64bits,  we need 256bits, resulting in 2.0bits per bit
 *   For 128bits, we need 192bits, resulting in 1.5bits per bit
 *   For 256bits, we need 320bits, resulting in 1.25bits per bit
 *
 */
template <size_t bits_ct, bool Align=true>
struct L0_NBitvector {
    std::vector<uint64_t>                      l0{0};
    std::vector<AlignedBitset<bits_ct, Align>> bits{{}};
    size_t totalLength{};


    L0_NBitvector() = default;
    L0_NBitvector(L0_NBitvector const&) = default;
    L0_NBitvector(L0_NBitvector&&) noexcept = default;

    // constructor accepting view to bools or already compact uint64_t
    template <std::ranges::sized_range range_t>
    L0_NBitvector(range_t&& _range) {
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
        l0.resize(totalLength/bits_ct + 1);
        bits.resize(totalLength/bits_ct + 1);
    }

    // the actual constructor, already receiving premade std::bitsets<N>
    template <std::ranges::sized_range range_t>
        requires std::same_as<std::ranges::range_value_t<range_t>, std::bitset<bits_ct>>
    L0_NBitvector(range_t&& _range) {
        auto _length = _range.size()*bits_ct;
        l0.resize(_length/bits_ct + 1);
        bits.resize(_length/bits_ct + 1);

        for (auto const& b : _range) {
            auto l0_id = totalLength / bits_ct;
            bits[l0_id].bits = b;
            totalLength += bits_ct;
            l0[l0_id+1] = l0[l0_id] + bits[l0_id].count();
        }
    }

    auto operator=(L0_NBitvector const&) -> L0_NBitvector& = default;
    auto operator=(L0_NBitvector&&) noexcept -> L0_NBitvector& = default;

    void reserve(size_t _length) {
        l0.reserve(_length/bits_ct + 1);
        bits.reserve(_length/bits_ct + 1);
    }

    void push_back(bool _value) {
        auto bitId         = totalLength % bits_ct;
        bits.back()[bitId] = _value;

        totalLength += 1;
        if (totalLength % bits_ct == 0) { // new l0-block
            l0.emplace_back(l0.back() + bits.back().count());
            bits.emplace_back();
        }
    }

    size_t size() const noexcept {
        return totalLength;
    }

    bool symbol(size_t idx) const noexcept {
        assert(idx < totalLength);
        auto bitId   = idx % bits_ct;
        auto blockId = idx / bits_ct;
        auto bit     = bits[blockId][bitId];
        return bit;
    }

    uint64_t rank(size_t idx) const noexcept {
        assert(idx <= totalLength);
        auto bitId = idx % bits_ct;
        auto l0Id  = idx / bits_ct;
        auto count = lshift_and_count(bits[l0Id].bits, bits_ct - bitId);
        auto r = l0[l0Id] + count;
        assert(r <= totalLength);
        return r;
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bits, l0, totalLength);
    }
};

using L0_64Bitvector  = L0_NBitvector<64>;
using L0_128Bitvector = L0_NBitvector<128>;
using L0_256Bitvector = L0_NBitvector<256>;
using L0_512Bitvector = L0_NBitvector<512>;
using L0_1024Bitvector = L0_NBitvector<1024>;
using L0_2048Bitvector = L0_NBitvector<2048>;

using L0_64BitvectorUA  = L0_NBitvector<64, false>;
using L0_128BitvectorUA = L0_NBitvector<128, false>;
using L0_256BitvectorUA = L0_NBitvector<256, false>;
using L0_512BitvectorUA = L0_NBitvector<512, false>;
using L0_1024BitvectorUA = L0_NBitvector<1024, false>;
using L0_2048BitvectorUA = L0_NBitvector<2048, false>;

static_assert(Bitvector_c<L0_64Bitvector>);
static_assert(Bitvector_c<L0_128Bitvector>);
static_assert(Bitvector_c<L0_256Bitvector>);
static_assert(Bitvector_c<L0_512Bitvector>);
static_assert(Bitvector_c<L0_1024Bitvector>);
static_assert(Bitvector_c<L0_2048Bitvector>);
static_assert(Bitvector_c<L0_64BitvectorUA>);
static_assert(Bitvector_c<L0_128BitvectorUA>);
static_assert(Bitvector_c<L0_256BitvectorUA>);
static_assert(Bitvector_c<L0_512BitvectorUA>);
static_assert(Bitvector_c<L0_1024BitvectorUA>);
static_assert(Bitvector_c<L0_2048BitvectorUA>);

}
