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
 * PairedL0L1_NBitvector a bit vector with only bits and blocks
 *
 */
template <size_t l1_bits_ct, size_t l0_bits_ct, bool Align=true, bool ShiftAndCount=false>
struct PairedL0L1_NBitvector {
    static_assert(l1_bits_ct < l0_bits_ct, "first level must be smaller than second level");
    static_assert(l0_bits_ct-l1_bits_ct <= std::numeric_limits<uint16_t>::max(), "l0_bits_ct can only hold up to uint16_t bits");
    std::vector<uint64_t> l0{0};
    std::vector<uint16_t> l1{0};
    std::vector<AlignedBitset<l1_bits_ct, Align>> bits{{}};
    size_t totalLength{};

    PairedL0L1_NBitvector() = default;
    PairedL0L1_NBitvector(PairedL0L1_NBitvector const&) = default;
    PairedL0L1_NBitvector(PairedL0L1_NBitvector&&) noexcept = default;

    // constructor accepting view to bools or already compact uint64_t
    template <std::ranges::sized_range range_t>
    PairedL0L1_NBitvector(range_t&& _range)
        : PairedL0L1_NBitvector{convertToBitsetView<l1_bits_ct>(std::forward<range_t>(_range))}
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
        l0.resize((totalLength+l0_bits_ct)/(l0_bits_ct*2) + 1);
        l1.resize((totalLength+l1_bits_ct)/(l1_bits_ct*2) + 1);
        bits.resize(totalLength/l1_bits_ct + 1);
    }

    // the actual constructor, already receiving premade std::bitsets<N>
    template <std::ranges::sized_range range_t>
        requires std::same_as<std::ranges::range_value_t<range_t>, std::bitset<l1_bits_ct>>
    PairedL0L1_NBitvector(range_t&& _range) {
        auto _length = _range.size()*l1_bits_ct;

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

    auto operator=(PairedL0L1_NBitvector const&) -> PairedL0L1_NBitvector& = default;
    auto operator=(PairedL0L1_NBitvector&&) noexcept -> PairedL0L1_NBitvector& = default;

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
    void serialize(Archive& ar) {
        ar(l0, l1, totalLength, bits);
    }
};
using PairedL0L1_64_4kBitvector   = PairedL0L1_NBitvector<64, 4096>;
using PairedL0L1_128_4kBitvector  = PairedL0L1_NBitvector<128, 4096>;
using PairedL0L1_256_4kBitvector  = PairedL0L1_NBitvector<256, 4096>;
using PairedL0L1_512_4kBitvector  = PairedL0L1_NBitvector<512, 4096>;
using PairedL0L1_1024_4kBitvector = PairedL0L1_NBitvector<1024, 4096>;
using PairedL0L1_2048_4kBitvector = PairedL0L1_NBitvector<2048, 4096>;

static_assert(BitVector_c<PairedL0L1_64_4kBitvector>);
static_assert(BitVector_c<PairedL0L1_128_4kBitvector>);
static_assert(BitVector_c<PairedL0L1_256_4kBitvector>);
static_assert(BitVector_c<PairedL0L1_512_4kBitvector>);
static_assert(BitVector_c<PairedL0L1_1024_4kBitvector>);
static_assert(BitVector_c<PairedL0L1_2048_4kBitvector>);

using PairedL0L1_64_64kBitvector   = PairedL0L1_NBitvector<64, 65536>;
using PairedL0L1_128_64kBitvector  = PairedL0L1_NBitvector<128, 65536>;
using PairedL0L1_256_64kBitvector  = PairedL0L1_NBitvector<256, 65536>;
using PairedL0L1_512_64kBitvector  = PairedL0L1_NBitvector<512, 65536>;
using PairedL0L1_1024_64kBitvector = PairedL0L1_NBitvector<1024, 65536>;
using PairedL0L1_2048_64kBitvector = PairedL0L1_NBitvector<2048, 65536>;

static_assert(BitVector_c<PairedL0L1_64_64kBitvector>);
static_assert(BitVector_c<PairedL0L1_128_64kBitvector>);
static_assert(BitVector_c<PairedL0L1_256_64kBitvector>);
static_assert(BitVector_c<PairedL0L1_512_64kBitvector>);
static_assert(BitVector_c<PairedL0L1_1024_64kBitvector>);
static_assert(BitVector_c<PairedL0L1_2048_64kBitvector>);

using PairedL0L1_64_64kBitvector_ShiftAndCount   = PairedL0L1_NBitvector<64, 65536, false, true>;
using PairedL0L1_128_64kBitvector_ShiftAndCount  = PairedL0L1_NBitvector<128, 65536, false, true>;
using PairedL0L1_256_64kBitvector_ShiftAndCount  = PairedL0L1_NBitvector<256, 65536, false, true>;
using PairedL0L1_512_64kBitvector_ShiftAndCount  = PairedL0L1_NBitvector<512, 65536, false, true>;
using PairedL0L1_1024_64kBitvector_ShiftAndCount = PairedL0L1_NBitvector<1024, 65536, false, true>;
using PairedL0L1_2048_64kBitvector_ShiftAndCount = PairedL0L1_NBitvector<2048, 65536, false, true>;

static_assert(BitVector_c<PairedL0L1_64_64kBitvector_ShiftAndCount>);
static_assert(BitVector_c<PairedL0L1_128_64kBitvector_ShiftAndCount>);
static_assert(BitVector_c<PairedL0L1_256_64kBitvector_ShiftAndCount>);
static_assert(BitVector_c<PairedL0L1_512_64kBitvector_ShiftAndCount>);
static_assert(BitVector_c<PairedL0L1_1024_64kBitvector_ShiftAndCount>);
static_assert(BitVector_c<PairedL0L1_2048_64kBitvector_ShiftAndCount>);


using PairedL0L1_64_64kBitvectorUA   = PairedL0L1_NBitvector<64, 65536, false>;
using PairedL0L1_128_64kBitvectorUA  = PairedL0L1_NBitvector<128, 65536, false>;
using PairedL0L1_256_64kBitvectorUA  = PairedL0L1_NBitvector<256, 65536, false>;
using PairedL0L1_512_64kBitvectorUA  = PairedL0L1_NBitvector<512, 65536, false>;
using PairedL0L1_1024_64kBitvectorUA = PairedL0L1_NBitvector<1024, 65536, false>;
using PairedL0L1_2048_64kBitvectorUA = PairedL0L1_NBitvector<2048, 65536, false>;

static_assert(BitVector_c<PairedL0L1_64_64kBitvectorUA>);
static_assert(BitVector_c<PairedL0L1_128_64kBitvectorUA>);
static_assert(BitVector_c<PairedL0L1_256_64kBitvectorUA>);
static_assert(BitVector_c<PairedL0L1_512_64kBitvectorUA>);
static_assert(BitVector_c<PairedL0L1_1024_64kBitvectorUA>);
static_assert(BitVector_c<PairedL0L1_2048_64kBitvectorUA>);


}
