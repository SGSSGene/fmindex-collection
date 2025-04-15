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
#include <limits>
#include <ranges>
#include <span>
#include <vector>

namespace fmindex_collection::bitvector {

/**
 * L0L1_NBitvector a bit vector with only bits and blocks
 *
 */
template <size_t l1_bits_ct, size_t l0_bits_ct, bool shift_and_count=false, bool Align=true>
struct L0L1_NBitvector {
    static_assert(l1_bits_ct < l0_bits_ct, "first level must be smaller than second level");
    static_assert(l0_bits_ct-l1_bits_ct <= std::numeric_limits<uint16_t>::max(), "l0_bits_ct can only hold up to uint16_t bits");
    std::vector<uint64_t> l0{0};
    std::vector<uint16_t> l1{0};
    std::vector<AlignedBitset<l1_bits_ct, Align>> bits{{}};
    size_t totalLength{};

    L0L1_NBitvector() = default;
    L0L1_NBitvector(L0L1_NBitvector const&) = default;
    L0L1_NBitvector(L0L1_NBitvector&&) noexcept = default;

    // constructor accepting view to bools or already compact uint64_t
    // converts it to a view that produces `std::bitset<l1_bits_ct>`
    template <std::ranges::sized_range range_t>
    L0L1_NBitvector(range_t&& _range)
        : L0L1_NBitvector{convertToBitsetView<l1_bits_ct>(std::forward<range_t>(_range))}
    {
        if constexpr (std::same_as<std::ranges::range_value_t<range_t>, uint64_t>) {
            totalLength = _range.size()*64;
        } else if constexpr (std::convertible_to<std::ranges::range_value_t<range_t>, bool>) {
            totalLength = _range.size();
        } else {
            []<bool b=false>() {
                static_assert(b, "Must be an uint64_t or convertible to bool");
            }();
        }
        l0.resize(totalLength/l0_bits_ct + 1);
        l1.resize(totalLength/l1_bits_ct + 1);
        bits.resize(totalLength/l1_bits_ct + 1);
    }


    // the actual constructor, already receiving premade std::bitsets<N>
    template <std::ranges::sized_range range_t>
        requires std::same_as<std::ranges::range_value_t<range_t>, std::bitset<l1_bits_ct>>
    L0L1_NBitvector(range_t&& _range) {
        auto _length = _range.size()*l1_bits_ct;
        l0.resize(_length/l0_bits_ct + 1);
        l1.resize(_length/l1_bits_ct + 1);
        bits.resize(_length/l1_bits_ct + 1);

        for (auto const& b : _range) {
            auto l1_id = totalLength / l1_bits_ct;
            auto l0_id = totalLength / l0_bits_ct;
            bits[l1_id].bits = b;
            totalLength += l1_bits_ct;
            l1[l1_id+1] = l1[l1_id] + bits[l1_id].count();
            if (totalLength % l0_bits_ct == 0) {
                l0[l0_id+1] = l0[l0_id] + l1[l1_id+1];
                l1[l1_id+1] = 0;
            }
        }
    }

    auto operator=(L0L1_NBitvector const&) -> L0L1_NBitvector& = default;
    auto operator=(L0L1_NBitvector&&) noexcept -> L0L1_NBitvector& = default;

    void reserve(size_t _length) {
        l0.reserve(_length/l0_bits_ct + 1);
        l1.reserve(_length/l1_bits_ct + 1);
        bits.reserve(_length/l1_bits_ct + 1);
    }

    void push_back(bool _value) {
        auto bitId         = totalLength % l1_bits_ct;
        bits.back()[bitId] = _value;

        totalLength += 1;
        if (totalLength % l1_bits_ct == 0) { // new l1-block
            l1.emplace_back(l1.back() + bits.back().count());
            bits.emplace_back();
            if (totalLength % l0_bits_ct == 0) { // new l0-block
                l0.emplace_back(l0.back() + l1.back());
                l1.back() = 0;
            }
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
        assert(idx <= totalLength);
        auto bitId = idx % (l1_bits_ct);
        auto l1Id = idx / l1_bits_ct;
        auto l0Id = idx / l0_bits_ct;
        assert(l1Id < bits.size());
        assert(l0Id < l0.size());

        auto count = [&]() {
            if constexpr (shift_and_count) {
                return (bits[l1Id].bits << (l1_bits_ct - bitId)).count();
            } else {
                return skip_first_or_last_n_bits_and_count(bits[l1Id].bits, bitId + l1_bits_ct);
            }
        }();

        auto r = l0[l0Id] + l1[l1Id] + count;
        assert(r <= idx);
        return r;
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(l0, l1, totalLength, bits);
    }
};
using L0L1_64_4kBitvector   = L0L1_NBitvector<64, 4096>;
using L0L1_128_4kBitvector  = L0L1_NBitvector<128, 4096>;
using L0L1_256_4kBitvector  = L0L1_NBitvector<256, 4096>;
using L0L1_512_4kBitvector  = L0L1_NBitvector<512, 4096>;
using L0L1_1024_4kBitvector = L0L1_NBitvector<1024, 4096>;
using L0L1_2048_4kBitvector = L0L1_NBitvector<2048, 4096>;

static_assert(BitVector_c<L0L1_64_4kBitvector>);
static_assert(BitVector_c<L0L1_128_4kBitvector>);
static_assert(BitVector_c<L0L1_256_4kBitvector>);
static_assert(BitVector_c<L0L1_512_4kBitvector>);
static_assert(BitVector_c<L0L1_1024_4kBitvector>);
static_assert(BitVector_c<L0L1_2048_4kBitvector>);

using L0L1_64_64kBitvector   = L0L1_NBitvector<64, 65536>;
using L0L1_128_64kBitvector  = L0L1_NBitvector<128, 65536>;
using L0L1_256_64kBitvector  = L0L1_NBitvector<256, 65536>;
using L0L1_512_64kBitvector  = L0L1_NBitvector<512, 65536>;
using L0L1_1024_64kBitvector = L0L1_NBitvector<1024, 65536>;
using L0L1_2048_64kBitvector = L0L1_NBitvector<2048, 65536>;

static_assert(BitVector_c<L0L1_64_64kBitvector>);
static_assert(BitVector_c<L0L1_128_64kBitvector>);
static_assert(BitVector_c<L0L1_256_64kBitvector>);
static_assert(BitVector_c<L0L1_512_64kBitvector>);
static_assert(BitVector_c<L0L1_1024_64kBitvector>);
static_assert(BitVector_c<L0L1_2048_64kBitvector>);

using L0L1_64_64kBitvector_ShiftAndCount   = L0L1_NBitvector<64, 65536, true>;
using L0L1_512_64kBitvector_ShiftAndCount  = L0L1_NBitvector<512, 65536, true>;
static_assert(BitVector_c<L0L1_64_64kBitvector_ShiftAndCount>);
static_assert(BitVector_c<L0L1_512_64kBitvector_ShiftAndCount>);

using L0L1_64_64kBitvectorUA   = L0L1_NBitvector<64, 65536, false, false>;
using L0L1_128_64kBitvectorUA  = L0L1_NBitvector<128, 65536, false, false>;
using L0L1_256_64kBitvectorUA  = L0L1_NBitvector<256, 65536, false, false>;
using L0L1_512_64kBitvectorUA  = L0L1_NBitvector<512, 65536, false, false>;
using L0L1_1024_64kBitvectorUA = L0L1_NBitvector<1024, 65536, false, false>;
using L0L1_2048_64kBitvectorUA = L0L1_NBitvector<2048, 65536, false, false>;

static_assert(BitVector_c<L0L1_64_64kBitvectorUA>);
static_assert(BitVector_c<L0L1_128_64kBitvectorUA>);
static_assert(BitVector_c<L0L1_256_64kBitvectorUA>);
static_assert(BitVector_c<L0L1_512_64kBitvectorUA>);
static_assert(BitVector_c<L0L1_1024_64kBitvectorUA>);
static_assert(BitVector_c<L0L1_2048_64kBitvectorUA>);

}
