// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../bitset_popcount.h"
#include "concepts.h"

#include <array>
#include <bitset>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <ranges>
#include <span>
#include <vector>

#if __has_include(<cereal/types/bitset.hpp>)
#include <cereal/types/bitset.hpp>
#endif

namespace fmindex_collection::bitvector {

/**
 * L1L2_NBitvector a bit vector with only bits and blocks
 *
 */
template <size_t l1_bits_ct, size_t l0_bits_ct, bool shift_and_count=false>
struct L1L2_NBitvector {
    static_assert(l1_bits_ct < l0_bits_ct, "first level must be smaller than second level");
    static_assert(l0_bits_ct-l1_bits_ct <= std::numeric_limits<uint16_t>::max(), "l0_bits_ct can only hold up to uint16_t bits");
    std::vector<uint64_t> l0{0};
    std::vector<uint16_t> l1{0};
    std::vector<std::bitset<l1_bits_ct>> bits{0};
    size_t totalLength{};
    bool finalized{};

    L1L2_NBitvector() = default;
    L1L2_NBitvector(L1L2_NBitvector const&) = default;
    L1L2_NBitvector(L1L2_NBitvector&&) noexcept = default;

    template <typename CB>
    L1L2_NBitvector(size_t length, CB cb)
        : L1L2_NBitvector{std::views::iota(size_t{}, length) | std::views::transform([&](size_t i) {
            return cb(i);
        })}
    {}

    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint8_t>
    L1L2_NBitvector(range_t&& _range) {
        reserve(_range.size());

        auto iter = _range.begin();
        while(totalLength < _range.size()) {
            push_back(*(iter++));
        }
    }

    auto operator=(L1L2_NBitvector const&) -> L1L2_NBitvector& = default;
    auto operator=(L1L2_NBitvector&&) noexcept -> L1L2_NBitvector& = default;

    void reserve(size_t _length) {
        l0.reserve(_length/(l0_bits_ct) + 2);
        l1.reserve(_length/(l1_bits_ct) + 2);
        bits.reserve(_length/l1_bits_ct + 2);
    }

    void finalize() const {
        const_cast<L1L2_NBitvector*>(this)->impl_finalize();
    }

    void impl_finalize() {
        if (finalized) return;
        finalized = true;

        size_t l0BlockCt = (totalLength + (l0_bits_ct)) / (l0_bits_ct);
        size_t l1BlockCt = l0BlockCt * (l0_bits_ct / l1_bits_ct);
        size_t inbitsCt  = l0BlockCt * ((l0_bits_ct) / l1_bits_ct);

        l0.resize(l0BlockCt);
        l1.resize(l1BlockCt);
        bits.resize(inbitsCt);

        constexpr size_t b1 = l1_bits_ct;
        constexpr size_t b0 = l0_bits_ct;

        constexpr size_t l1_block_ct = b0 / b1;


        size_t l0_acc{};
        // walk through all superblocks
        for (size_t l0I{0}; l0I < l0BlockCt; ++l0I) {
            // update l0
            l0[l0I] = l0_acc;

            // right part
            size_t acc = 0;
            for (size_t i{0}; i < l1_block_ct; ++i) {
                auto idx = l0I*l1_block_ct + i;
                l1[idx] = acc;
                acc += bits[l0I*l1_block_ct + i].count();
            }
            l0_acc += acc;
        }
    }

    void push_back(bool _value) {
        if (finalized) throw std::runtime_error{"DS is finalized, can not call push_back after calling `rank`"};
        if (_value) {
            auto bitId         = totalLength % l1_bits_ct;
            bits.back()[bitId] = _value;
        }
        totalLength += 1;
        if (totalLength % l1_bits_ct == 0) { // next bit will require a new in-bits block
            bits.emplace_back();
        }
    }

    size_t size() const noexcept {
        return totalLength;
    }

    bool symbol(size_t idx) const noexcept {
        assert(idx < totalLength);
        auto bitId = idx % l1_bits_ct;
        auto l1Id  = idx / l1_bits_ct;
        assert(l1Id < bits.size());
        auto bit = bits[l1Id][bitId];
        return bit;
    }

    uint64_t rank(size_t idx) const noexcept {
        finalize();
        assert(idx <= totalLength);
        auto bitId = idx % (l1_bits_ct);
        auto l1Id = idx / l1_bits_ct;
        auto l0Id = idx / l0_bits_ct;
        assert(l1Id < bits.size());
        assert(l0Id < l0.size());

        auto count = [&]() {
            if constexpr (shift_and_count) {
                return (bits[l1Id] << (l1_bits_ct - bitId)).count();
            } else {
                return skip_first_or_last_n_bits_and_count(bits[l1Id], bitId + l1_bits_ct);
            }
        }();

        auto r = l0[l0Id] + l1[l1Id] + count;
        assert(r <= idx);
        return r;
    }

    template <typename Archive>
    void save(Archive& ar) const {
        finalize();
        ar(l0, l1, totalLength);
        saveBV(bits, ar);
    }

    template <typename Archive>
    void load(Archive& ar) {
        ar(l0, l1, totalLength);
        loadBV(bits, ar);
    }

};
using L1L2_64_4kBitvector   = L1L2_NBitvector<64, 4096>;
using L1L2_128_4kBitvector  = L1L2_NBitvector<128, 4096>;
using L1L2_256_4kBitvector  = L1L2_NBitvector<256, 4096>;
using L1L2_512_4kBitvector  = L1L2_NBitvector<512, 4096>;
using L1L2_1024_4kBitvector = L1L2_NBitvector<1024, 4096>;
using L1L2_2048_4kBitvector = L1L2_NBitvector<2048, 4096>;

static_assert(BitVector_c<L1L2_64_4kBitvector>);
static_assert(BitVector_c<L1L2_128_4kBitvector>);
static_assert(BitVector_c<L1L2_256_4kBitvector>);
static_assert(BitVector_c<L1L2_512_4kBitvector>);
static_assert(BitVector_c<L1L2_1024_4kBitvector>);
static_assert(BitVector_c<L1L2_2048_4kBitvector>);

using L1L2_64_64kBitvector   = L1L2_NBitvector<64, 65536>;
using L1L2_128_64kBitvector  = L1L2_NBitvector<128, 65536>;
using L1L2_256_64kBitvector  = L1L2_NBitvector<256, 65536>;
using L1L2_512_64kBitvector  = L1L2_NBitvector<512, 65536>;
using L1L2_1024_64kBitvector = L1L2_NBitvector<1024, 65536>;
using L1L2_2048_64kBitvector = L1L2_NBitvector<2048, 65536>;

static_assert(BitVector_c<L1L2_64_64kBitvector>);
static_assert(BitVector_c<L1L2_128_64kBitvector>);
static_assert(BitVector_c<L1L2_256_64kBitvector>);
static_assert(BitVector_c<L1L2_512_64kBitvector>);
static_assert(BitVector_c<L1L2_1024_64kBitvector>);
static_assert(BitVector_c<L1L2_2048_64kBitvector>);

using L1L2_64_64kBitvector_ShiftAndCount   = L1L2_NBitvector<64, 65536, true>;
using L1L2_512_64kBitvector_ShiftAndCount  = L1L2_NBitvector<512, 65536, true>;
static_assert(BitVector_c<L1L2_64_64kBitvector_ShiftAndCount>);
static_assert(BitVector_c<L1L2_512_64kBitvector_ShiftAndCount>);

}
