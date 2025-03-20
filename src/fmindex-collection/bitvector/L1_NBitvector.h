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
    std::vector<uint64_t> superblocks{0};
    std::vector<std::bitset<bits_ct>> bits{0};
    size_t totalLength{};

    static std::array<std::bitset<bits_ct>, bits_ct+1> const masks;/* = []() {
        auto m = std::array<std::bitset<bits_ct>, bits_ct+1>{};
        auto a = std::bitset<bits_ct>{};

        for (size_t i{0}; i < bits_ct; ++i) {
            a.set(i);
            m[i+1] = a;
        }
        return m;
    }();*/

    L1_NBitvector() = default;
    L1_NBitvector(L1_NBitvector const&) = default;
    L1_NBitvector(L1_NBitvector&&) noexcept = default;

    template <typename CB>
    L1_NBitvector(size_t length, CB cb)
        : L1_NBitvector{std::views::iota(size_t{}, length) | std::views::transform([&](size_t i) {
            return cb(i);
        })}
    {}

    template <typename CB>
    struct Layer {
        CB get;

        auto operator[](size_t i) -> decltype(auto) {
            return get(i);
        }
    };

    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint8_t>
    L1_NBitvector(range_t&& _range) {
        reserve(_range.size());

        auto _length = _range.size();
        superblocks.resize(_length/bits_ct + 1);
        bits.resize(_length/bits_ct + 1);

        auto l0 = Layer{[&](size_t i) -> uint64_t& {
            return superblocks[i];
        }};
        auto l1 = Layer{[&](size_t i) -> decltype(auto) {
            return bits[i];
        }};

        for (auto iter = _range.begin(); iter != _range.end(); ++iter) {
            // run only if full block
            auto restBits = std::min(size_t{bits_ct}, _range.size() - totalLength);

            // concatenate next full block
            auto l0_id = totalLength / bits_ct;
            auto& bits = l1[l0_id];
            bits = bool(*iter);
            if (restBits == bits_ct) {
                for (size_t j{0}; j < bits_ct/64; ++j) {
                    uint64_t v{};
                    for (size_t i{(j==0)?1ull:0ull}; i < 64; ++i) {
                        bool value = *(++iter);
                        v |= (uint64_t{value} << i);
                    }
                    bits = bits | (std::bitset<bits_ct>{v} << (j*64));
                }
            } else {
                for (size_t i{1}; i < restBits; ++i) {
                    bool value = *(++iter);
                    bits = bits | (std::bitset<bits_ct>{value} << i);
                }
            }

            totalLength += restBits;

            // abort - if block not full,
            if (restBits < bits_ct) {
                break;
            }
            l0[l0_id+1] = l0[l0_id] + bits.count();
        }
    }

    auto operator=(L1_NBitvector const&) -> L1_NBitvector& = default;
    auto operator=(L1_NBitvector&&) noexcept -> L1_NBitvector& = default;

    void reserve(size_t _length) {
        superblocks.reserve(_length/bits_ct + 1);
        bits.reserve(_length/bits_ct + 1);
    }

    void push_back(bool _value) {
        auto bitId         = totalLength % bits_ct;
        bits.back()[bitId] = _value;

        totalLength += 1;
        if (totalLength % bits_ct == 0) { // new superblock
            superblocks.emplace_back(superblocks.back() + bits.back().count());
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

        auto bitcount = (bits[superblockId] & masks[bitId]).count();
//        auto maskedBits = (bits[superblockId] << (bits_ct - bitId));
//        auto bitcount   = maskedBits.count();

        return superblocks[superblockId]
                + bitcount;
    }

    template <typename Archive>
    void save(Archive& ar) const {
        ar(superblocks, totalLength);
        size_t ct = bits.size();
        ar(ct);

        auto mask = std::bitset<bits_ct>{~uint64_t{0}};

        for (auto const& b : bits) {
            // saving in 64bit blocks
            for (size_t i{0}; i < bits_ct; i += 64) {
                auto v = ((b >> i) & mask).to_ullong();
                ar(v);
            }
        }
        //ar(superblocks, bits, totalLength);
    }

    template <typename Archive>
    void load(Archive& ar) {
        ar(superblocks, totalLength);
        size_t ct{};
        ar(ct);

        bits.resize(ct);

        for (size_t j{0}; j < ct; ++j) {
            auto b = std::bitset<bits_ct>{};
            // saving in 64bit blocks
            for (size_t i{0}; i < bits_ct; i += 64) {
                uint64_t v{};
                ar(v);
                b = b | (std::bitset<bits_ct>{v} << i);
            }
            bits[j] = b;
        }
        //ar(superblocks, bits, totalLength);
    }
};

template <size_t bits_ct>
std::array<std::bitset<bits_ct>, bits_ct+1> const L1_NBitvector<bits_ct>::masks = []() {
    auto m = std::array<std::bitset<bits_ct>, bits_ct+1>{};
    auto a = std::bitset<bits_ct>{};

    for (size_t i{0}; i < bits_ct; ++i) {
        a.set(i);
        m[i+1] = a;
    }
    return m;
}();

using L1_64Bitvector  = L1_NBitvector<64>;
using L1_128Bitvector = L1_NBitvector<128>;
using L1_256Bitvector = L1_NBitvector<256>;
using L1_512Bitvector = L1_NBitvector<512>;
using L1_1024Bitvector = L1_NBitvector<1024>;
using L1_2048Bitvector = L1_NBitvector<2048>;

static_assert(BitVector_c<L1_64Bitvector>);
static_assert(BitVector_c<L1_128Bitvector>);
static_assert(BitVector_c<L1_256Bitvector>);
static_assert(BitVector_c<L1_512Bitvector>);
static_assert(BitVector_c<L1_1024Bitvector>);
static_assert(BitVector_c<L1_2048Bitvector>);

}
