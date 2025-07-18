// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../FixedSuccinctVector2.h"
#include "../bitset_popcount.h"
#include "concepts.h"

#include <array>
#include <bit>
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

namespace fmc::bitvector {

#if !__clang__ || __clang_major__ >= 19 // !TODO workaround (weird optimization bug?)

/**
 * CompactPairedL1L2_NBitvector a bit vector with only bits and blocks
 *
 */
template <size_t l1_bits_ct, size_t l0_bits_ct>
struct CompactPairedL1L2_NBitvector {
    static_assert(l1_bits_ct < l0_bits_ct, "first level must be smaller than second level");
    static constexpr size_t word_width_l0 = 64;
    static constexpr size_t word_width_l1 = std::bit_width(l0_bits_ct-l1_bits_ct);
    std::vector<uint64_t> l0{0};
    FixedSuccinctVector2<word_width_l1> l1{0};
    std::vector<std::bitset<l1_bits_ct>> bits{0};
    size_t totalLength{};
    bool finalized{};

    CompactPairedL1L2_NBitvector() = default;
    CompactPairedL1L2_NBitvector(CompactPairedL1L2_NBitvector const&) = default;
    CompactPairedL1L2_NBitvector(CompactPairedL1L2_NBitvector&&) noexcept = default;

    template <typename CB>
    CompactPairedL1L2_NBitvector(size_t length, CB cb)
        : CompactPairedL1L2_NBitvector{std::views::iota(size_t{}, length) | std::views::transform([&](size_t i) {
            return cb(i);
        })}
    {}

    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint8_t>
    CompactPairedL1L2_NBitvector(range_t&& _range) {
        reserve(_range.size());

        auto iter = _range.begin();
        while(totalLength < _range.size()) {
            push_back(*(iter++));
        }
    }

    auto operator=(CompactPairedL1L2_NBitvector const&) -> CompactPairedL1L2_NBitvector& = default;
    auto operator=(CompactPairedL1L2_NBitvector&&) noexcept -> CompactPairedL1L2_NBitvector& = default;

    void reserve(size_t _length) {
        l0.reserve(_length/(l1_bits_ct*2) + 2);
        l1.reserve(_length/(l0_bits_ct*2) + 2);
        bits.reserve(_length/l1_bits_ct + 2);
    }

    void finalize() const {
        const_cast<CompactPairedL1L2_NBitvector*>(this)->impl_finalize();
    }

    void impl_finalize() {
        if (finalized) return;
        finalized = true;

        size_t l0BlockCt = (totalLength + (l0_bits_ct*2)) / (l0_bits_ct*2);
        size_t l1BlockCt = l0BlockCt * (l0_bits_ct / l1_bits_ct);
        size_t inbitsCt  = l0BlockCt * ((l0_bits_ct*2) / l1_bits_ct);

        l0.resize(l0BlockCt);
        l1.resize(l1BlockCt);
        bits.resize(inbitsCt);

        constexpr size_t b1 = l1_bits_ct;
        constexpr size_t b0 = l0_bits_ct;

        constexpr size_t l1_block_ct = b0 / b1;


        size_t l0_acc{};
        // walk through all superblocks
        for (size_t l0I{0}; l0I < l0BlockCt; ++l0I) {
            // walk left to right and set l1 values (as if they are the begining of a superblock)

            // left part
            size_t acc{};
            for (size_t i{0}; i < l1_block_ct; ++i) {
                acc += bits[l0I*l1_block_ct*2 + i].count();
                if (i % 2 == 0) {
                    l1.set(l0I*l1_block_ct + i/2, acc);
                }
            }
            l0_acc += acc;
            // update l0 (reached center)
            l0[l0I] = l0_acc;
//            l0.set(l0I, l0_acc);//acc + [&]() { if(l0I > 0) return l0[l0I-1]; return size_t{0};}();
            // walk backwards through left part and revert l0
            for (size_t i{0}; i < l1_block_ct; ++i) {
                if (i % 2 == 0) {
                    auto idx = l0I*l1_block_ct + i/2;
                    l1.set(idx, acc - l1.at(idx));
                }
            }

            // right part
            acc = 0;
            for (size_t i{l1_block_ct}; i < l1_block_ct*2; ++i) {
                acc += bits[l0I*l1_block_ct*2 + i].count();
                if (i % 2 == 0) {
                    auto idx = l0I*l1_block_ct + i/2;
                    l1.set(idx, acc);
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
        assert(l0Id/2 < l0.size());

        int64_t right_l1 = (l1Id%2)*2-1;
        int64_t right_l0 = (l0Id%2)*2-1;

        auto count = skip_first_or_last_n_bits_and_count(bits[l1Id], bitId);

        auto r = l0[l0Id/2] + right_l0 * l1[l1Id/2] + right_l1 * count;
        assert(r <= idx);
        return r;
    }

    template <typename Archive>
    void save(Archive& ar) const {
        finalize();
        ar(l0, l1, totalLength, finalized);
        saveBV(bits, ar);
    }

    template <typename Archive>
    void load(Archive& ar) {
        ar(l0, l1, totalLength, finalized);
        loadBV(bits, ar);
    }
};

using CompactPairedL1L2_64_4kBitvector   = CompactPairedL1L2_NBitvector<64, 4096>;
using CompactPairedL1L2_128_4kBitvector  = CompactPairedL1L2_NBitvector<128, 4096>;
using CompactPairedL1L2_256_4kBitvector  = CompactPairedL1L2_NBitvector<256, 4096>;
using CompactPairedL1L2_512_4kBitvector  = CompactPairedL1L2_NBitvector<512, 4096>;
using CompactPairedL1L2_1024_4kBitvector = CompactPairedL1L2_NBitvector<1024, 4096>;
using CompactPairedL1L2_2048_4kBitvector = CompactPairedL1L2_NBitvector<2048, 4096>;

static_assert(Bitvector_c<CompactPairedL1L2_64_4kBitvector>);
static_assert(Bitvector_c<CompactPairedL1L2_128_4kBitvector>);
static_assert(Bitvector_c<CompactPairedL1L2_256_4kBitvector>);
static_assert(Bitvector_c<CompactPairedL1L2_512_4kBitvector>);
static_assert(Bitvector_c<CompactPairedL1L2_1024_4kBitvector>);
static_assert(Bitvector_c<CompactPairedL1L2_2048_4kBitvector>);

using CompactPairedL1L2_64_64kBitvector   = CompactPairedL1L2_NBitvector<64, 65536>;
using CompactPairedL1L2_128_64kBitvector  = CompactPairedL1L2_NBitvector<128, 65536>;
using CompactPairedL1L2_256_64kBitvector  = CompactPairedL1L2_NBitvector<256, 65536>;
using CompactPairedL1L2_512_64kBitvector  = CompactPairedL1L2_NBitvector<512, 65536>;
using CompactPairedL1L2_1024_64kBitvector = CompactPairedL1L2_NBitvector<1024, 65536>;
using CompactPairedL1L2_2048_64kBitvector = CompactPairedL1L2_NBitvector<2048, 65536>;

static_assert(Bitvector_c<CompactPairedL1L2_64_64kBitvector>);
static_assert(Bitvector_c<CompactPairedL1L2_128_64kBitvector>);
static_assert(Bitvector_c<CompactPairedL1L2_256_64kBitvector>);
static_assert(Bitvector_c<CompactPairedL1L2_512_64kBitvector>);
static_assert(Bitvector_c<CompactPairedL1L2_1024_64kBitvector>);
static_assert(Bitvector_c<CompactPairedL1L2_2048_64kBitvector>);
#endif

}
