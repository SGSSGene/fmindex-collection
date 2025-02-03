// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../FixedSuccinctVector.h"
#include "../bitset_popcount.h"
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
 * CompactDoubleL1L2_NBitvector a bit vector with only bits and blocks
 *
 */
template <size_t bits_ct, size_t bits_ct2>
struct CompactDoubleL1L2_NBitvector {
    static constexpr size_t word_width_l0 = 64;
    static constexpr size_t word_width_l1 = std::bit_width(bits_ct2);
    FixedSuccinctVector<word_width_l0> l0{0, 0};
    FixedSuccinctVector<word_width_l1> l1{0};
//    std::vector<uint64_t> l0{0, 0};
//    std::vector<uint16_t> l1{0};
    std::vector<std::bitset<bits_ct>> bits{0};
    size_t totalLength{};

    CompactDoubleL1L2_NBitvector() = default;
    CompactDoubleL1L2_NBitvector(CompactDoubleL1L2_NBitvector const&) = default;
    CompactDoubleL1L2_NBitvector(CompactDoubleL1L2_NBitvector&&) noexcept = default;

    template <typename CB>
    CompactDoubleL1L2_NBitvector(size_t length, CB cb)
        : CompactDoubleL1L2_NBitvector{std::views::iota(size_t{}, length) | std::views::transform([&](size_t i) {
            return cb(i);
        })}
    {}

    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint8_t>
    CompactDoubleL1L2_NBitvector(range_t&& _range) {
        reserve(_range.size());

        auto iter = _range.begin();

/*        size_t const loop64  = _range.size() / (bits_ct*2);

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
        totalLength = bits_ct * (loop64*2);*/
        size_t acc{};
        while(totalLength < _range.size()) {
            auto v = *(iter++);
            acc += v;
            push_back(v);
            if (rank(totalLength) != acc) {
                auto s = rank(totalLength);
                (void)s;
            }
            assert(rank(totalLength) == acc);
        }
    }

    auto operator=(CompactDoubleL1L2_NBitvector const&) -> CompactDoubleL1L2_NBitvector& = default;
    auto operator=(CompactDoubleL1L2_NBitvector&&) noexcept -> CompactDoubleL1L2_NBitvector& = default;

    void reserve(size_t _length) {
        l0.reserve((_length+1)/(bits_ct) + 1);
        l1.reserve((_length+1)/(bits_ct2*2) + 1);
        bits.reserve(_length/bits_ct + 1);
    }

    void push_back(bool _value) {
        if (_value) {
            auto bitId         = totalLength % bits_ct;
            bits.back()[bitId] = _value;
            assert(l1.back() < 4096);
            { auto t = l0.back(); l0.pop_back(); l0.push_back(t + _value); }
            { auto t = l1.back(); l1.pop_back(); l1.push_back(t + _value); }
        }
        totalLength += 1;
        if (totalLength % bits_ct2 == 0) { // new l0 block
            l0.emplace_back(l0.back());
            l1.pop_back();
            l1.push_back(0);
            assert(l1.back() == 0);
            bits.emplace_back();
        } else if (totalLength % (bits_ct*2) == bits_ct) { // new l1 block
            l1.emplace_back(l1.back());
            bits.emplace_back();
        } else if (totalLength % (bits_ct*2) == 0) { // new in-bits block
            bits.emplace_back();
        }
    }

    size_t size() const noexcept {
        return totalLength;
    }

    bool symbol(size_t idx) const noexcept {
        assert(idx < totalLength);
        auto bitId = idx % bits_ct;
        auto l1Id  = idx / bits_ct;
        auto bit = bits[l1Id][bitId];
        return bit;
    }

    uint64_t rank(size_t idx) const noexcept {
        assert(idx <= totalLength);
        auto bitId = idx % (bits_ct*2);
        auto l1Id = idx / bits_ct;
        auto l0Id = idx / bits_ct2;

        auto right_l1 = (l1Id%2);

        auto count = signed_rshift_and_count(bits[l1Id], bitId);

        auto r = l0[l0Id] + l1[l1Id/2] + (right_l1*2-1) * count;
        assert(r <= idx);
        return r;
    }

    template <typename Archive>
    void save(Archive& ar) const {
        ar(l0, l1, totalLength);
        saveBV(bits, ar);
    }

    template <typename Archive>
    void load(Archive& ar) {
        ar(l0, l1, totalLength);
        loadBV(bits, ar);
    }

};
using CompactDoubleL1L2_64_4kBitvector   = CompactDoubleL1L2_NBitvector<64, 4096>;
using CompactDoubleL1L2_128_4kBitvector  = CompactDoubleL1L2_NBitvector<128, 4096>;
using CompactDoubleL1L2_256_4kBitvector  = CompactDoubleL1L2_NBitvector<256, 4096>;
using CompactDoubleL1L2_512_4kBitvector  = CompactDoubleL1L2_NBitvector<512, 4096>;
using CompactDoubleL1L2_1024_4kBitvector = CompactDoubleL1L2_NBitvector<1024, 4096>;
using CompactDoubleL1L2_2048_4kBitvector = CompactDoubleL1L2_NBitvector<2048, 4096>;

static_assert(BitVector_c<CompactDoubleL1L2_64_4kBitvector>);
static_assert(BitVector_c<CompactDoubleL1L2_128_4kBitvector>);
static_assert(BitVector_c<CompactDoubleL1L2_256_4kBitvector>);
static_assert(BitVector_c<CompactDoubleL1L2_512_4kBitvector>);
static_assert(BitVector_c<CompactDoubleL1L2_1024_4kBitvector>);
static_assert(BitVector_c<CompactDoubleL1L2_2048_4kBitvector>);

using CompactDoubleL1L2_64_64kBitvector   = CompactDoubleL1L2_NBitvector<64, 65536>;
using CompactDoubleL1L2_128_64kBitvector  = CompactDoubleL1L2_NBitvector<128, 65536>;
using CompactDoubleL1L2_256_64kBitvector  = CompactDoubleL1L2_NBitvector<256, 65536>;
using CompactDoubleL1L2_512_64kBitvector  = CompactDoubleL1L2_NBitvector<512, 65536>;
using CompactDoubleL1L2_1024_64kBitvector = CompactDoubleL1L2_NBitvector<1024, 65536>;
using CompactDoubleL1L2_2048_64kBitvector = CompactDoubleL1L2_NBitvector<2048, 65536>;

static_assert(BitVector_c<CompactDoubleL1L2_64_64kBitvector>);
static_assert(BitVector_c<CompactDoubleL1L2_128_64kBitvector>);
static_assert(BitVector_c<CompactDoubleL1L2_256_64kBitvector>);
static_assert(BitVector_c<CompactDoubleL1L2_512_64kBitvector>);
static_assert(BitVector_c<CompactDoubleL1L2_1024_64kBitvector>);
static_assert(BitVector_c<CompactDoubleL1L2_2048_64kBitvector>);


}
