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
 * InterleavedPairedL0L1_NBitvector a bit vector with only bits and blocks
 *
 */
template <size_t l1_bits_ct, size_t l0_bits_ct, bool Align=true>
struct InterleavedPairedL0L1_NBitvector {
    static_assert(l1_bits_ct < l0_bits_ct, "first level must be smaller than second level");
    static_assert(l0_bits_ct-l1_bits_ct <= std::numeric_limits<uint16_t>::max(), "l0_bits_ct can only hold up to uint16_t bits");

    struct Entries {
        uint64_t l0;
        std::array<uint16_t, l0_bits_ct / l1_bits_ct> l1;

        template <typename Archive>
        void serialize(Archive& ar) {
            ar(l0, l1);
        }
    };
    std::vector<Entries> l0l1{Entries{{0}, {0}}};
    std::vector<AlignedBitset<l1_bits_ct, Align>> bits{{}};
    size_t totalLength{};
    bool finalized{};

    InterleavedPairedL0L1_NBitvector() = default;
    InterleavedPairedL0L1_NBitvector(InterleavedPairedL0L1_NBitvector const&) = default;
    InterleavedPairedL0L1_NBitvector(InterleavedPairedL0L1_NBitvector&&) noexcept = default;

    template <typename CB>
    InterleavedPairedL0L1_NBitvector(size_t length, CB cb)
        : InterleavedPairedL0L1_NBitvector{std::views::iota(size_t{}, length) | std::views::transform([&](size_t i) {
            return cb(i);
        })}
    {}

    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint8_t>
    InterleavedPairedL0L1_NBitvector(range_t&& _range) {
        reserve(_range.size());

        auto iter = _range.begin();
        while(totalLength < _range.size()) {
            push_back(*(iter++));
        }
    }

    auto operator=(InterleavedPairedL0L1_NBitvector const&) -> InterleavedPairedL0L1_NBitvector& = default;
    auto operator=(InterleavedPairedL0L1_NBitvector&&) noexcept -> InterleavedPairedL0L1_NBitvector& = default;

    void reserve(size_t _length) {
        l0l1.reserve(_length/(l0_bits_ct*2) + 2);
        bits.reserve(_length/l1_bits_ct + 2);
    }

    void finalize() const {
        const_cast<InterleavedPairedL0L1_NBitvector*>(this)->impl_finalize();
    }

    void impl_finalize() {
        if (finalized) return;
        finalized = true;

        size_t l0BlockCt = (totalLength + (l0_bits_ct*2)) / (l0_bits_ct*2);
        //size_t l1BlockCt = l0BlockCt * (l0_bits_ct / l1_bits_ct);
        size_t inbitsCt  = l0BlockCt * ((l0_bits_ct*2) / l1_bits_ct);

        l0l1.resize(l0BlockCt);
        bits.resize(inbitsCt);

        constexpr size_t b1 = l1_bits_ct;
        constexpr size_t b0 = l0_bits_ct;

        constexpr size_t l1_block_ct = b0 / b1;
//        static_assert(l1_block_ct == 128, "must be equal to the number of element in the Entries struct");


        size_t l0_acc{};
        // walk through all superblocks
        for (size_t l0I{0}; l0I < l0BlockCt; ++l0I) {
            // walk left to right and set l1 values (as if they are the begining of a superblock)

            // left part
            size_t acc{};
            for (size_t i{0}; i < l1_block_ct; ++i) {
                acc += bits[l0I*l1_block_ct*2 + i].count();
                if (i % 2 == 0) {
                    l0l1[l0I].l1[i/2] = acc;
//                    l1[l0I*l1_block_ct + i/2] = acc;
                }
            }
            l0_acc += acc;
            // update l0 (reached center)
            l0l1[l0I].l0 = l0_acc;
//            l0[l0I] = l0_acc;//acc + [&]() { if(l0I > 0) return l0[l0I-1]; return size_t{0};}();
            // walk backwards through left part and revert l0
            for (size_t i{0}; i < l1_block_ct; ++i) {
                if (i % 2 == 0) {
//                    auto idx = l0I*l1_block_ct + i/2;
                    l0l1[l0I].l1[i/2] = acc - l0l1[l0I].l1[i/2];
//                    l1[idx] = acc - l1[idx];
                }
            }

            // right part
            acc = 0;
            for (size_t i{l1_block_ct}; i < l1_block_ct*2; ++i) {
                acc += bits[l0I*l1_block_ct*2 + i].count();
                if (i % 2 == 0) {
  //                  auto idx = l0I*l1_block_ct + i/2;
                    l0l1[l0I].l1[i/2] = acc;
//                    l1[idx] = acc;
                }
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
        assert(l1Id/2 < bits.size());
        auto bit = bits[l1Id][bitId];
        return bit;
    }

    uint64_t rank(size_t idx) const noexcept {
        finalize();
        assert(idx <= totalLength);
        auto bitId = idx % (l1_bits_ct*2);
        auto l1Id = idx / l1_bits_ct;
        auto l0Id = idx / l0_bits_ct;
        assert(l1Id/2 < bits.size());
        assert(l0Id/2 < l0l1.size());

        int64_t right_l1 = (l1Id%2)*2-1;
        int64_t right_l0 = (l0Id%2)*2-1;

        auto count = skip_first_or_last_n_bits_and_count(bits[l1Id].bits, bitId);
        auto& l0 = l0l1[l0Id/2];
        auto l1  = l0.l1[(idx - l0Id * l0_bits_ct) / l1_bits_ct / 2];
        auto r = l0.l0 + right_l0 * l1 + right_l1 * count;
//        auto r = l0[l0Id/2] + right_l0 * l1[l1Id/2] + right_l1 * count;
        assert(r <= idx);
        return r;
    }

    template <typename Archive>
    void save(Archive& ar) const {
        finalize();
        ar(l0l1, totalLength, bits);
    }

    template <typename Archive>
    void load(Archive& ar) {
        ar(l0l1, totalLength, bits);
    }

};

using InterleavedPairedL0L1_64_4kBitvector   = InterleavedPairedL0L1_NBitvector<64, 4096>;
using InterleavedPairedL0L1_128_4kBitvector  = InterleavedPairedL0L1_NBitvector<128, 4096>;
using InterleavedPairedL0L1_256_4kBitvector  = InterleavedPairedL0L1_NBitvector<256, 4096>;
using InterleavedPairedL0L1_512_4kBitvector  = InterleavedPairedL0L1_NBitvector<512, 4096>;
using InterleavedPairedL0L1_1024_4kBitvector = InterleavedPairedL0L1_NBitvector<1024, 4096>;
using InterleavedPairedL0L1_2048_4kBitvector = InterleavedPairedL0L1_NBitvector<2048, 4096>;

static_assert(Bitvector_c<InterleavedPairedL0L1_64_4kBitvector>);
static_assert(Bitvector_c<InterleavedPairedL0L1_128_4kBitvector>);
static_assert(Bitvector_c<InterleavedPairedL0L1_256_4kBitvector>);
static_assert(Bitvector_c<InterleavedPairedL0L1_512_4kBitvector>);
static_assert(Bitvector_c<InterleavedPairedL0L1_1024_4kBitvector>);
static_assert(Bitvector_c<InterleavedPairedL0L1_2048_4kBitvector>);

using InterleavedPairedL0L1_64_64kBitvector   = InterleavedPairedL0L1_NBitvector<64, 65536>;
using InterleavedPairedL0L1_128_64kBitvector  = InterleavedPairedL0L1_NBitvector<128, 65536>;
using InterleavedPairedL0L1_256_64kBitvector  = InterleavedPairedL0L1_NBitvector<256, 65536>;
using InterleavedPairedL0L1_512_64kBitvector  = InterleavedPairedL0L1_NBitvector<512, 65536>;
using InterleavedPairedL0L1_1024_64kBitvector = InterleavedPairedL0L1_NBitvector<1024, 65536>;
using InterleavedPairedL0L1_2048_64kBitvector = InterleavedPairedL0L1_NBitvector<2048, 65536>;

static_assert(Bitvector_c<InterleavedPairedL0L1_64_64kBitvector>);
static_assert(Bitvector_c<InterleavedPairedL0L1_128_64kBitvector>);
static_assert(Bitvector_c<InterleavedPairedL0L1_256_64kBitvector>);
static_assert(Bitvector_c<InterleavedPairedL0L1_512_64kBitvector>);
static_assert(Bitvector_c<InterleavedPairedL0L1_1024_64kBitvector>);
static_assert(Bitvector_c<InterleavedPairedL0L1_2048_64kBitvector>);

using InterleavedPairedL0L1_64_64kBitvectorUA   = InterleavedPairedL0L1_NBitvector<64, 65536, false>;
using InterleavedPairedL0L1_128_64kBitvectorUA  = InterleavedPairedL0L1_NBitvector<128, 65536, false>;
using InterleavedPairedL0L1_256_64kBitvectorUA  = InterleavedPairedL0L1_NBitvector<256, 65536, false>;
using InterleavedPairedL0L1_512_64kBitvectorUA  = InterleavedPairedL0L1_NBitvector<512, 65536, false>;
using InterleavedPairedL0L1_1024_64kBitvectorUA = InterleavedPairedL0L1_NBitvector<1024, 65536, false>;
using InterleavedPairedL0L1_2048_64kBitvectorUA = InterleavedPairedL0L1_NBitvector<2048, 65536, false>;

static_assert(Bitvector_c<InterleavedPairedL0L1_64_64kBitvectorUA>);
static_assert(Bitvector_c<InterleavedPairedL0L1_128_64kBitvectorUA>);
static_assert(Bitvector_c<InterleavedPairedL0L1_256_64kBitvectorUA>);
static_assert(Bitvector_c<InterleavedPairedL0L1_512_64kBitvectorUA>);
static_assert(Bitvector_c<InterleavedPairedL0L1_1024_64kBitvectorUA>);
static_assert(Bitvector_c<InterleavedPairedL0L1_2048_64kBitvectorUA>);

}
