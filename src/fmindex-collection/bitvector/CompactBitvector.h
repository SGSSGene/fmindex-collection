// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "concepts.h"

#include <bitset>
#include <cassert>
#include <ranges>
#include <vector>

namespace fmindex_collection::bitvector {

/**
 * CompactBitvector with interleaved superblocks, blocks and bits
 *
 * - Each group consist of 384bits, divided into 6 blocks.
 * - Each block uses 9bits to represents a value (6*9bits = 54bits).
 *   The last 10bits are padding bits, not used for any thing.
 * - Superblock consist of a single 64bit number
 *
 *   For 384bits, we need 512bits, or 1.333bits to save a single bit
 */
struct CompactBitvector {
    struct alignas(64) Superblock {
        uint64_t superBlockEntry{};
        uint64_t blockEntries{};
        std::array<uint64_t, 6> bits{};

        uint64_t rank(size_t idx) const noexcept {
            assert(idx < 384);

            auto blockId = idx >> 6;
            auto block = uint64_t{0b111111111} & (blockEntries >> (blockId * 9));
            auto keep = (uint64_t{idx} & 63);
            if (keep == 0) return superBlockEntry + block;

            auto maskedBits = bits[blockId] << (64 - keep);
            auto bitcount   = std::bitset<64>{maskedBits}.count();

            auto total = superBlockEntry + block + bitcount;
            return total;
        }

        bool symbol(size_t idx) const noexcept {
            assert(idx < 384);

            auto blockId = idx >> 6;
            auto bitId = idx & 63;
            return bits[blockId] & (uint64_t{1} << bitId);
        }

        void setBlock(size_t blockId, size_t value) {
            blockEntries = blockEntries & ~(uint64_t{0b111111111} << blockId*9);
            blockEntries = blockEntries | (uint64_t{value} << blockId*9);
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar(superBlockEntry, blockEntries, bits);
        }
    };

    static constexpr size_t Sigma = 2;

    std::vector<Superblock> superblocks{Superblock{}};
    size_t                  totalLength{};

    template <typename CB>
    CompactBitvector(size_t length, CB cb)
        : CompactBitvector{std::views::iota(size_t{}, length) | std::views::transform([&](size_t i) {
            return cb(i);
        })}
    {}

    //!TODO helper structures, when building the vector
    uint64_t sblock_acc{};
    uint64_t block_acc{};

    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint8_t>
    CompactBitvector(range_t&& _range) {

        reserve(_range.size());

        auto iter = _range.begin();
        while (totalLength < _range.size()) {
            push_back(*(iter++));
        }
    }

    CompactBitvector() = default;
    CompactBitvector(CompactBitvector const&) = default;
    CompactBitvector(CompactBitvector&&) noexcept = default;
    auto operator=(CompactBitvector const&) -> CompactBitvector& = default;
    auto operator=(CompactBitvector&&) noexcept -> CompactBitvector& = default;

    void reserve(size_t _length) {
        superblocks.reserve((_length+1)/(64*6) + 1);
    }

    void push_back(bool _value) {
        if (_value) {
            auto blockId      = (totalLength >>  6) % 6;
            auto bitId        = totalLength &  63;

            auto& bits = superblocks.back().bits[blockId];
            bits = bits | (1ull << bitId);

            block_acc  += 1;
            sblock_acc += 1;
        }
        totalLength += 1;
        if (totalLength % 64 == 0) { // new block
            if (totalLength % 384 == 0) { // new super block + new block
                superblocks.emplace_back();
                superblocks.back().superBlockEntry = sblock_acc;
                block_acc = 0;
            } else {
                superblocks.back().setBlock((totalLength % 384) / 64, block_acc);
            }
        }
    }


    size_t size() const noexcept {
        return totalLength;
    }

    bool symbol(size_t idx) const noexcept {
        auto superblockId = idx / 384;
        auto bitId        = idx % 384;
        return superblocks[superblockId].symbol(bitId);
    }

    uint64_t rank(size_t idx) const noexcept {
        assert(idx <= size());
        auto superblockId = idx / 384;
        auto bitId        = idx % 384;
        auto v = superblocks[superblockId].rank(bitId);
        assert(v <= idx);
        return v;
    }



    template <typename Archive>
    void serialize(Archive& ar) {
        ar(totalLength, superblocks, sblock_acc, block_acc);
    }
};
static_assert(BitVector_c<CompactBitvector>);

}
