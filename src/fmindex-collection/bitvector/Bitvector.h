// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "concepts.h"
#include "../bitset_popcount.h"

#include <array>
#include <bit>
#include <bitset>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <ranges>
#include <span>
#include <vector>

namespace fmc::bitvector {

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
    std::vector<uint64_t> superblocks{0};
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

    template <typename CB>
    struct Layer {
        CB get;

        auto operator[](size_t i) -> auto& {
            return get(i);
        }
    };

    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint8_t>
    Bitvector(range_t&& _range) {
        reserve(_range.size());

        auto _length = _range.size();
        superblocks.resize(_length/256+1);
        blocks.resize(_length/64+1);
        bits.resize(_length/64+1);

        auto l0 = Layer{[&](size_t i) -> uint64_t& {
            return superblocks[i];
        }};
        auto l1 = Layer{[&](size_t i) -> uint8_t& {
            return blocks[i];
        }};
        auto l2 = Layer{[&](size_t i) -> uint64_t& {
            return bits[i];
        }};

        for (auto iter = _range.begin(); iter != _range.end(); ++iter) {
            // run only if full block
            auto restBits = std::min(size_t{64}, _range.size() - totalLength);

            // concatenate next full block
            uint64_t bits = *iter;
            for (size_t i{1}; i < restBits; ++i) {
                bool value = *(++iter);
                bits = bits | (uint64_t{value} << i);
            }

            // update bits and blocks
            auto l0_id = totalLength / 256;
            auto l1_id = totalLength / 64;
            l2[l1_id] = bits;

            totalLength += restBits;
            // abort - if block not full,
            if (restBits < 64) {
                break;
            }

            auto ct = size_t{l1[l1_id]} + std::popcount(bits);
            l1[l1_id+1] = ct;
            // check if next superblock is full
            if (totalLength % 256 == 0) {
                l0[l0_id+1] = l0[l0_id] + ct;
                l1[l1_id+1] = 0;
            }
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
        }
        totalLength += 1;
        if (totalLength % 64 == 0) { // new block
            auto ct = size_t{blocks.back()} + std::popcount(bits.back());
            blocks.emplace_back(ct);
            bits.emplace_back();
            if (totalLength % 256 == 0) { // new super block + new block
                superblocks.emplace_back(superblocks.back() + ct);
                blocks.back() = 0;
            }
        }
    }

    size_t size() const noexcept {
        return totalLength;
    }

    bool symbol(size_t idx) const noexcept {
        assert(idx <= totalLength);
        auto bitId        = idx % 64;
        auto blockId      = idx / 64;
        auto bit = (bits[blockId] >> bitId) & 1;
        return bit;
    }

    uint64_t rank(size_t idx) const noexcept {
        assert(idx <= totalLength);
        auto bitId        = idx % 64;
        auto blockId      = idx / 64;
        auto superblockId = blockId / 4;
        if (idx % 64 == 0) {
            auto r = superblocks[superblockId] + blocks[blockId];
            assert(r <= idx);
            return r;
        }

        auto maskedBits = (bits[blockId] << (64 - bitId));
        auto bitcount   = std::bitset<64>{maskedBits}.count();

        auto r = superblocks[superblockId]
                + blocks[blockId]
                + bitcount;
        assert(r <= idx);
        return r;
    }

    template <typename Archive>
    void serialize(this auto&& self, Archive& ar) {
        ar(self.superblocks, self.blocks, self.bits, self.totalLength);
    }

    static size_t estimateSize(size_t totalSize) {
        auto bits_for_l0_blocks = (totalSize/256) * 64;
        auto bits_for_l1_blocks = (totalSize/8) * 8;
        return totalSize + bits_for_l0_blocks + bits_for_l1_blocks;
    }

};

static_assert(Bitvector_c<Bitvector>);

}
