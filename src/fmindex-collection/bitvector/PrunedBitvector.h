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
#include <iostream>

namespace fmc::bitvector {

/**
 * Like Bitvector but collapses to a few bytes if everything is zero
 */
struct PrunedBitvector {
    std::vector<uint64_t> superblocks{0, 0};
    std::vector<uint8_t>  blocks{0};
    std::vector<uint64_t> bits{0};
    size_t totalLength{};

    bool allZeros{true};
    size_t totalLengthAllZeros{};

    PrunedBitvector() = default;
    PrunedBitvector(PrunedBitvector const&) = default;
    PrunedBitvector(PrunedBitvector&&) noexcept = default;

    template <typename CB>
    PrunedBitvector(size_t length, CB cb)
        : PrunedBitvector{std::views::iota(size_t{}, length) | std::views::transform([&](size_t i) {
            return cb(i);
        })}
    {}

    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint8_t>
    PrunedBitvector(range_t&& _range) {
//        reserve(_range.size());

        // check if all zeros
        allZeros = [&]() {
            for (auto v : _range) {
                if (v != 0) {
                    return false;
                }
            }
            return true;
        }();
        if (allZeros) {
            totalLengthAllZeros = _range.size();
            return;
        }

        auto iter = _range.begin();

        size_t const loop64  = _range.size() / 64;
        for (size_t l64{}; l64 < loop64; ++l64) {
            for (size_t i{}; i < 64; ++i) {
                bool value = *(iter++);
                if (value) {
                    auto bitId = i;
                    auto& tbits  = bits.back();
                    tbits        = tbits | (uint64_t{1} << bitId);
                    superblocks.back() += 1;
                }
            }
            totalLength += 64;
            blocks.emplace_back(superblocks.back());
            bits.emplace_back();
            if (totalLength % 256 == 0) {
                superblocks.back() += superblocks[superblocks.size()-2];
                superblocks.emplace_back();
                blocks.back() = 0;
            }
        }
        while(totalLength < _range.size()) {
            push_back(*(iter++));
        }
    }

    auto operator=(PrunedBitvector const&) -> PrunedBitvector& = default;
    auto operator=(PrunedBitvector&&) noexcept -> PrunedBitvector& = default;

    void reserve(size_t _length) {
        superblocks.reserve((_length+1)/(64*4) + 1);
        blocks.reserve((_length+1)/64 + 1);
        bits.reserve(_length/64 + 1);
    }

    void push_back(bool _value) {
        if (!allZeros) {
            push_backImpl(_value);
        } else if (_value) {
            for (size_t i{0}; i < totalLengthAllZeros; ++i) {
                push_backImpl(0);
            }
            push_backImpl(1);
            allZeros = false;
        } else {
            totalLengthAllZeros = totalLengthAllZeros + 1;
        }
    }

    void push_backImpl(bool _value) {
        if (_value) {
            auto bitId   = totalLength % 64;
            auto& tbits  = bits.back();
            tbits        = tbits | (uint64_t{1} << bitId);
            superblocks.back() += 1;
        }
        totalLength += 1;
        if (totalLength % 64 == 0) { // new block
            if (totalLength % 256 == 0) { // new super block + new block
                superblocks.back() += superblocks[superblocks.size()-2];
                superblocks.emplace_back();
                blocks.emplace_back();
                bits.emplace_back();
            } else {
                assert(superblocks.back() < 256);
                blocks.emplace_back(uint8_t(superblocks.back()));
                bits.emplace_back();
            }
        }
    }


    size_t size() const noexcept {
        if (allZeros) return totalLengthAllZeros;
        return totalLength;
    }

    bool symbol(size_t idx) const noexcept {
        if (allZeros) return 0;
        auto bitId        = idx % 64;
        auto blockId      = idx / 64;
        auto bit = (bits[blockId] >> bitId) & 1;
        return bit;
    }

    uint64_t rank(size_t idx) const noexcept {
        if (allZeros) return 0;
        auto bitId        = idx % 64;
        auto blockId      = idx / 64;
        auto superblockId = blockId / 4;
        if (idx % 64 == 0) return superblocks[superblockId] + blocks[blockId];

        auto maskedBits = (bits[blockId] << (64 - bitId));
        auto bitcount   = std::bitset<64>{maskedBits}.count();

        return superblocks[superblockId]
                + blocks[blockId]
                + bitcount;
    }

    template <typename Archive>
    void serialize(this auto&& self, Archive& ar) {
        ar(self.superblocks, self.blocks, self.bits, self.totalLength, self.allZeros, self.totalLengthAllZeros);
    }
};

static_assert(Bitvector_c<PrunedBitvector>);

}
