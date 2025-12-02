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
 * PairedBitvector2L a bit vector with only bits and blocks
 *
 */
template <size_t l1_bits_ct, size_t l0_bits_ct, bool Align=true, bool ShiftAndCount=false>
struct PairedBitvector2L {
    static_assert(l1_bits_ct < l0_bits_ct, "first level must be smaller than second level");
    static_assert(l0_bits_ct-l1_bits_ct <= std::numeric_limits<uint16_t>::max(), "l0_bits_ct can only hold up to uint16_t bits");
    std::vector<uint64_t> l0{0};
    std::vector<uint16_t> l1{0};
    std::vector<AlignedBitset<l1_bits_ct, Align>> bits{{}};
    size_t totalLength{};

    PairedBitvector2L() = default;
    PairedBitvector2L(PairedBitvector2L const&) = default;
    PairedBitvector2L(PairedBitvector2L&&) noexcept = default;

    // constructor accepting view to bools or already compact uint64_t
    template <std::ranges::sized_range range_t>
    PairedBitvector2L(range_t&& _range) {
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
        l0.resize((totalLength+l0_bits_ct)/(l0_bits_ct*2) + 1);
        l1.resize((totalLength+l1_bits_ct)/(l1_bits_ct*2) + 1);
        bits.resize(totalLength/l1_bits_ct + 1);
    }

    // the actual constructor, already receiving premade std::bitsets<N>
    template <std::ranges::sized_range range_t>
        requires std::same_as<std::ranges::range_value_t<range_t>, std::bitset<l1_bits_ct>>
    PairedBitvector2L(range_t&& _range) {
        auto _length = static_cast<size_t>(_range.size()*l1_bits_ct);

        l0.resize((_length+l0_bits_ct)/(l0_bits_ct*2) + 1);
        l1.resize((_length+l1_bits_ct)/(l1_bits_ct*2) + 1);
        bits.resize(_length/l1_bits_ct + 1);

        for (auto const& b : _range) {
            auto l0_id   = totalLength / (l0_bits_ct*2);
            auto l1_id   = totalLength / (l1_bits_ct*2);
            auto bits_id = totalLength / l1_bits_ct;
            bits[bits_id].bits = b;

            auto isLeftL0 = totalLength % (l0_bits_ct*2) < l0_bits_ct;
            auto isLeftL1 = totalLength % (l1_bits_ct*2) < l1_bits_ct;

            auto count = bits[bits_id].count();
            if (isLeftL0) {
                l0[l0_id] += count;
                auto startL1Id = l0_id * l0_bits_ct / l1_bits_ct;

                for (size_t i{startL1Id}; i < l1_id + !isLeftL1; ++i) {
                    l1[i] += count;
                }
            } else {
                l0[l0_id+1] += count;
                if (isLeftL1) {
                    l1[l1_id] += count;
                    l1[l1_id+1] = l1[l1_id];
                } else {
                    l1[l1_id+1] += count;
                }
            }
            totalLength += l1_bits_ct;

            // requires extension to next blocks?
            if (isLeftL0) {
                // switches from left to right in next block
                if (totalLength % (l0_bits_ct*2) >= l0_bits_ct) {
                    l0[l0_id+1] = l0[l0_id];
                }
            } else {
                // switches from right to left in next block
                if (totalLength % (l0_bits_ct*2) < l0_bits_ct) {
                    l1[l1_id+1] = 0;
                }
            }
        }
    }

    auto operator=(PairedBitvector2L const&) -> PairedBitvector2L& = default;
    auto operator=(PairedBitvector2L&&) noexcept -> PairedBitvector2L& = default;

    void reserve(size_t _length) {
        l0.reserve(_length/(l0_bits_ct*2) + 2);
        l1.reserve(_length/(l1_bits_ct*2) + 2);
        bits.reserve(_length/l1_bits_ct + 2);
    }

    void push_back(bool _value) {
        auto bitId         = totalLength % l1_bits_ct;
        bits.back()[bitId] = _value;

        auto l0_id   = totalLength / (l0_bits_ct*2);
        auto l1_id   = totalLength / (l1_bits_ct*2);

        size_t count = _value;

        auto isLeftL0 = totalLength % (l0_bits_ct*2) < l0_bits_ct;
        auto isLeftL1 = totalLength % (l1_bits_ct*2) < l1_bits_ct;

        // add to the appropiate counters
        if (isLeftL0) {
            l0.back() += count;

            auto startL1Id = l0_id * l0_bits_ct / l1_bits_ct;
            for (size_t i{startL1Id}; i < l1_id + !isLeftL1; ++i) {
                l1[i] += count;
            }
        } else {
            l0.back() += count;
            l1.back() += count;
        }

        totalLength += 1;
        // requires extension to next blocks?
        if (isLeftL0) {
            // switches l0 from left to right in next block
            if (totalLength % (l0_bits_ct*2) >= l0_bits_ct) {
                l0.emplace_back(l0.back());
            }
            if (totalLength % (l1_bits_ct*2) == l1_bits_ct) {
                l1.emplace_back();
            }
        } else {
            if (totalLength % (l1_bits_ct*2) == l1_bits_ct) {
                l1.emplace_back(l1.back());
            }
            // switches l0 from right to left in next block
            if (totalLength % (l0_bits_ct*2) < l0_bits_ct) {
                l1.back() = 0;
            }
        }
        if (totalLength % l1_bits_ct == 0) {
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
        assert(idx <= totalLength);
        auto bitId = idx % (l1_bits_ct*2);
        auto l1Id = idx / l1_bits_ct;
        auto l0Id = idx / l0_bits_ct;
        assert(l1Id < bits.size());
        assert(l1Id/2 < l1.size());
        assert(l0Id/2 < l0.size());

        int64_t right_l1 = (l1Id%2)*2-1;
        int64_t right_l0 = (l0Id%2)*2-1;

        auto count = [&]() -> size_t {
            if constexpr (ShiftAndCount) {
                if (bitId > l1_bits_ct) {
                    size_t i = l1_bits_ct - bitId;
                    return (bits[l1Id].bits << i).count();
                } else {
                    return (bits[l1Id].bits >> bitId).count();
                }
            } else {
                return skip_first_or_last_n_bits_and_count(bits[l1Id].bits, bitId);
            }
        }();

        auto r = l0[l0Id/2] + right_l0 * l1[l1Id/2] + right_l1 * count;
        assert(r <= idx);
        return r;
    }

    template <typename Archive>
    void serialize(this auto&& self, Archive& ar) {
        ar(self.l0, self.l1, self.totalLength, self.bits);
    }
};
using PairedBitvector2L_64_4k   = PairedBitvector2L<64, 4096>;
using PairedBitvector2L_128_4k  = PairedBitvector2L<128, 4096>;
using PairedBitvector2L_256_4k  = PairedBitvector2L<256, 4096>;
using PairedBitvector2L_512_4k  = PairedBitvector2L<512, 4096>;
using PairedBitvector2L_1024_4k = PairedBitvector2L<1024, 4096>;
using PairedBitvector2L_2048_4k = PairedBitvector2L<2048, 4096>;

static_assert(Bitvector_c<PairedBitvector2L_64_4k>);
static_assert(Bitvector_c<PairedBitvector2L_128_4k>);
static_assert(Bitvector_c<PairedBitvector2L_256_4k>);
static_assert(Bitvector_c<PairedBitvector2L_512_4k>);
static_assert(Bitvector_c<PairedBitvector2L_1024_4k>);
static_assert(Bitvector_c<PairedBitvector2L_2048_4k>);

using PairedBitvector2L_64_64k   = PairedBitvector2L<64, 65536>;
using PairedBitvector2L_128_64k  = PairedBitvector2L<128, 65536>;
using PairedBitvector2L_256_64k  = PairedBitvector2L<256, 65536>;
using PairedBitvector2L_512_64k  = PairedBitvector2L<512, 65536>;
using PairedBitvector2L_1024_64k = PairedBitvector2L<1024, 65536>;
using PairedBitvector2L_2048_64k = PairedBitvector2L<2048, 65536>;

static_assert(Bitvector_c<PairedBitvector2L_64_64k>);
static_assert(Bitvector_c<PairedBitvector2L_128_64k>);
static_assert(Bitvector_c<PairedBitvector2L_256_64k>);
static_assert(Bitvector_c<PairedBitvector2L_512_64k>);
static_assert(Bitvector_c<PairedBitvector2L_1024_64k>);
static_assert(Bitvector_c<PairedBitvector2L_2048_64k>);

using PairedBitvector2L_64_64k_ShiftAndCount   = PairedBitvector2L<64, 65536, false, true>;
using PairedBitvector2L_128_64k_ShiftAndCount  = PairedBitvector2L<128, 65536, false, true>;
using PairedBitvector2L_256_64k_ShiftAndCount  = PairedBitvector2L<256, 65536, false, true>;
using PairedBitvector2L_512_64k_ShiftAndCount  = PairedBitvector2L<512, 65536, false, true>;
using PairedBitvector2L_1024_64k_ShiftAndCount = PairedBitvector2L<1024, 65536, false, true>;
using PairedBitvector2L_2048_64k_ShiftAndCount = PairedBitvector2L<2048, 65536, false, true>;

static_assert(Bitvector_c<PairedBitvector2L_64_64k_ShiftAndCount>);
static_assert(Bitvector_c<PairedBitvector2L_128_64k_ShiftAndCount>);
static_assert(Bitvector_c<PairedBitvector2L_256_64k_ShiftAndCount>);
static_assert(Bitvector_c<PairedBitvector2L_512_64k_ShiftAndCount>);
static_assert(Bitvector_c<PairedBitvector2L_1024_64k_ShiftAndCount>);
static_assert(Bitvector_c<PairedBitvector2L_2048_64k_ShiftAndCount>);


using PairedBitvector2L_64_64kUA   = PairedBitvector2L<64, 65536, false>;
using PairedBitvector2L_128_64kUA  = PairedBitvector2L<128, 65536, false>;
using PairedBitvector2L_256_64kUA  = PairedBitvector2L<256, 65536, false>;
using PairedBitvector2L_512_64kUA  = PairedBitvector2L<512, 65536, false>;
using PairedBitvector2L_1024_64kUA = PairedBitvector2L<1024, 65536, false>;
using PairedBitvector2L_2048_64kUA = PairedBitvector2L<2048, 65536, false>;

static_assert(Bitvector_c<PairedBitvector2L_64_64kUA>);
static_assert(Bitvector_c<PairedBitvector2L_128_64kUA>);
static_assert(Bitvector_c<PairedBitvector2L_256_64kUA>);
static_assert(Bitvector_c<PairedBitvector2L_512_64kUA>);
static_assert(Bitvector_c<PairedBitvector2L_1024_64kUA>);
static_assert(Bitvector_c<PairedBitvector2L_2048_64kUA>);


}
