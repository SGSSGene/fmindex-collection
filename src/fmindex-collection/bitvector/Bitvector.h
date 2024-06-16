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

namespace fmindex_collection::bitvector {

/**
 * Bitvector with interleaved superblocks, blocks and bits
 *
 * - Each group consist of 256bits, divided into 4 blocks.
 * - Each block uses 8bits to represents a value (4*8bits = 32bits).
 * - Superblock consist of a single 64bit number
 *
 *   For 256bits, we need 352bits, or 1.375bits to save a single bit
 */
struct Bitvector {
    std::vector<uint64_t> superblocks{0, 0};
    std::vector<uint8_t>  blocks{0};
    std::vector<uint64_t> bits{0};
    size_t totalLength{};

    Bitvector() = default;
    Bitvector(Bitvector const&) = default;
    Bitvector(Bitvector&&) noexcept = default;

    template <typename CB>
    Bitvector(size_t length, CB cb)
        : Bitvector{std::views::iota(size_t{}, length) | std::views::transform([&](size_t i) {
            return cb(i);
        })}
    {}

    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint8_t>
    Bitvector(range_t&& _range) {
        reserve(_range.size());

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

    auto operator=(Bitvector const&) -> Bitvector& = default;
    auto operator=(Bitvector&&) noexcept -> Bitvector& = default;

    void reserve(size_t _length) {
        superblocks.reserve((_length+1)/(64*4) + 1);
        blocks.reserve((_length+1)/64 + 1);
        bits.reserve(_length/64 + 1);
    }

    void push_back(bool _value) {
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
        return totalLength;
    }

    bool symbol(size_t idx) const noexcept {
        auto bitId        = idx % 64;
        auto blockId      = idx / 64;
        auto bit = (bits[blockId] >> bitId) & 1;
        return bit;
    }

    uint64_t rank(size_t idx) const noexcept {
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
    void serialize(Archive& ar) {
        ar(superblocks, blocks, bits, totalLength);
    }
};

static_assert(BitVector_c<Bitvector>);

}
