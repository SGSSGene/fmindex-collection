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

namespace fmc::bitvector {

/**
 * Bitvector2L a bit vector with only bits and blocks
 *
 */
template <size_t l1_bits_ct, size_t l0_bits_ct, bool shift_and_count=false, bool Align=true>
struct Bitvector2L {
    static_assert(l1_bits_ct < l0_bits_ct, "first level must be smaller than second level");
    static_assert(l0_bits_ct-l1_bits_ct <= std::numeric_limits<uint16_t>::max(), "l0_bits_ct can only hold up to uint16_t bits");
    std::vector<uint64_t> l0{0};
    std::vector<uint16_t> l1{0};
    std::vector<AlignedBitset<l1_bits_ct, Align>> bits{{}};
    size_t totalLength{};

    Bitvector2L() = default;
    Bitvector2L(Bitvector2L const&) = default;
    Bitvector2L(Bitvector2L&&) noexcept = default;

    // constructor accepting view to bools or already compact uint64_t
    // converts it to a view that produces `std::bitset<l1_bits_ct>`
    template <std::ranges::sized_range range_t>
    Bitvector2L(range_t&& _range) {
        auto _size = _range.size();
        if constexpr (std::same_as<std::ranges::range_value_t<range_t>, uint64_t>) {
            *this = {std::forward<range_t>(_range) | view_as_bitset<l1_bits_ct>};
            totalLength = _size*64;
        } else if constexpr (std::convertible_to<std::ranges::range_value_t<range_t>, bool>) {
            *this = {std::forward<range_t>(_range) | view_bool_as_uint64 | view_as_bitset<l1_bits_ct>};
            totalLength = _size;
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
    Bitvector2L(range_t&& _range) {
        auto _length = _range.size()*l1_bits_ct;
        l0.resize(_length/l0_bits_ct + 1);
        l1.resize(_length/l1_bits_ct + 1);
        bits.resize(_length/l1_bits_ct + 1);

        uint64_t l1_a{};
        for (auto const& b : _range) {
            auto l1_id = totalLength / l1_bits_ct;
            auto l0_id = totalLength / l0_bits_ct;
            bits[l1_id].bits = b;
            totalLength += l1_bits_ct;
            l1_a = l1_a + bits[l1_id].count();
            l1[l1_id+1] = l1_a;
            if (totalLength % l0_bits_ct == 0) {
                l0[l0_id+1] = l0[l0_id] + l1_a;
                l1[l1_id+1] = 0;
                l1_a = 0;
            }
        }
    }

    auto operator=(Bitvector2L const&) -> Bitvector2L& = default;
    auto operator=(Bitvector2L&&) noexcept -> Bitvector2L& = default;

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

    uint64_t gotoMarkingFwd(size_t idx) const {
        assert(idx < totalLength);
        while (!symbol(idx)) {
            idx += 1;
        }
        return idx;
    }

    uint64_t gotoMarkingBwd(size_t idx) const {
        assert(idx < totalLength);
        while (!symbol(idx)) {
            idx -= 1;
        }
        return idx;
    }

    template <typename Archive>
    void serialize(this auto&& self, Archive& ar) {
        ar(self.l0, self.l1, self.totalLength, self.bits);
    }

    static size_t estimateSize(size_t totalSize) {
        auto bits_for_l0_blocks = (totalSize/l0_bits_ct + 1) * 64;
        auto bits_for_l1_blocks = (totalSize/l1_bits_ct + 1) * 16;
        auto bits = (totalSize/l1_bits_ct+1) * l1_bits_ct;
        return bits + bits_for_l0_blocks + bits_for_l1_blocks;
    }
};

using Bitvector2L_64_4k   = Bitvector2L<64, 4096>;
using Bitvector2L_128_4k  = Bitvector2L<128, 4096>;
using Bitvector2L_256_4k  = Bitvector2L<256, 4096>;
using Bitvector2L_512_4k  = Bitvector2L<512, 4096>;
using Bitvector2L_1024_4k = Bitvector2L<1024, 4096>;
using Bitvector2L_2048_4k = Bitvector2L<2048, 4096>;

static_assert(Bitvector_c<Bitvector2L_64_4k>);
static_assert(Bitvector_c<Bitvector2L_128_4k>);
static_assert(Bitvector_c<Bitvector2L_256_4k>);
static_assert(Bitvector_c<Bitvector2L_512_4k>);
static_assert(Bitvector_c<Bitvector2L_1024_4k>);
static_assert(Bitvector_c<Bitvector2L_2048_4k>);

using Bitvector2L_64_64k   = Bitvector2L<64, 65536>;
using Bitvector2L_128_64k  = Bitvector2L<128, 65536>;
using Bitvector2L_256_64k  = Bitvector2L<256, 65536>;
using Bitvector2L_512_64k  = Bitvector2L<512, 65536>;
using Bitvector2L_1024_64k = Bitvector2L<1024, 65536>;
using Bitvector2L_2048_64k = Bitvector2L<2048, 65536>;

static_assert(Bitvector_c<Bitvector2L_64_64k>);
static_assert(Bitvector_c<Bitvector2L_128_64k>);
static_assert(Bitvector_c<Bitvector2L_256_64k>);
static_assert(Bitvector_c<Bitvector2L_512_64k>);
static_assert(Bitvector_c<Bitvector2L_1024_64k>);
static_assert(Bitvector_c<Bitvector2L_2048_64k>);

using Bitvector2L_64_64k_ShiftAndCount   = Bitvector2L<64, 65536, true>;
using Bitvector2L_512_64k_ShiftAndCount  = Bitvector2L<512, 65536, true>;
static_assert(Bitvector_c<Bitvector2L_64_64k_ShiftAndCount>);
static_assert(Bitvector_c<Bitvector2L_512_64k_ShiftAndCount>);

using Bitvector2L_64_64kUA   = Bitvector2L<64, 65536, false, false>;
using Bitvector2L_128_64kUA  = Bitvector2L<128, 65536, false, false>;
using Bitvector2L_256_64kUA  = Bitvector2L<256, 65536, false, false>;
using Bitvector2L_512_64kUA  = Bitvector2L<512, 65536, false, false>;
using Bitvector2L_1024_64kUA = Bitvector2L<1024, 65536, false, false>;
using Bitvector2L_2048_64kUA = Bitvector2L<2048, 65536, false, false>;

static_assert(Bitvector_c<Bitvector2L_64_64kUA>);
static_assert(Bitvector_c<Bitvector2L_128_64kUA>);
static_assert(Bitvector_c<Bitvector2L_256_64kUA>);
static_assert(Bitvector_c<Bitvector2L_512_64kUA>);
static_assert(Bitvector_c<Bitvector2L_1024_64kUA>);
static_assert(Bitvector_c<Bitvector2L_2048_64kUA>);

}
