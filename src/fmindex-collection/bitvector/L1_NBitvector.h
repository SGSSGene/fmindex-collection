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

namespace fmindex_collection::bitvector {

/**
 * L1_NBitvector a bit vector with only bits and blocks
 *
 *   For 64bits,  we need 256bits, resulting in 2.0bits per bit
 *   For 128bits, we need 192bits, resulting in 1.5bits per bit
 *   For 256bits, we need 320bits, resulting in 1.25bits per bit
 *
 */
template <size_t bits_ct, bool Align=true>
struct L1_NBitvector {
    std::vector<uint64_t>                      superblocks{0};
    std::vector<AlignedBitset<bits_ct, Align>> bits{{}};
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
            auto& bits = l1[l0_id].bits;
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
        auto count = lshift_and_count(bits[superblockId].bits, bits_ct - bitId);
        return superblocks[superblockId]
                + count;
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(superblocks, totalLength, bits);
    }
};

using L1_64Bitvector  = L1_NBitvector<64>;
using L1_128Bitvector = L1_NBitvector<128>;
using L1_256Bitvector = L1_NBitvector<256>;
using L1_512Bitvector = L1_NBitvector<512>;
using L1_1024Bitvector = L1_NBitvector<1024>;
using L1_2048Bitvector = L1_NBitvector<2048>;

using L1_64BitvectorUA  = L1_NBitvector<64, false>;
using L1_128BitvectorUA = L1_NBitvector<128, false>;
using L1_256BitvectorUA = L1_NBitvector<256, false>;
using L1_512BitvectorUA = L1_NBitvector<512, false>;
using L1_1024BitvectorUA = L1_NBitvector<1024, false>;
using L1_2048BitvectorUA = L1_NBitvector<2048, false>;

static_assert(BitVector_c<L1_64Bitvector>);
static_assert(BitVector_c<L1_128Bitvector>);
static_assert(BitVector_c<L1_256Bitvector>);
static_assert(BitVector_c<L1_512Bitvector>);
static_assert(BitVector_c<L1_1024Bitvector>);
static_assert(BitVector_c<L1_2048Bitvector>);
static_assert(BitVector_c<L1_64BitvectorUA>);
static_assert(BitVector_c<L1_128BitvectorUA>);
static_assert(BitVector_c<L1_256BitvectorUA>);
static_assert(BitVector_c<L1_512BitvectorUA>);
static_assert(BitVector_c<L1_1024BitvectorUA>);
static_assert(BitVector_c<L1_2048BitvectorUA>);

}
