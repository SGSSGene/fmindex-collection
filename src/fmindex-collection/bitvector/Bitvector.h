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
#include <span>
#include <vector>

namespace fmindex_collection {
namespace bitvector {

/**
 * Bitvector with interleaved superblocks, blocks and bits
 *
 * - Each group consist of 384bits, divided into 6 blocks.
 * - Each block uses 9bits to represents a value (6*9bits = 54bits).
 *   The last 10bits are padding bits, not used for any thing.
 * - Superblock consist of a single 64bit number
 *
 *   For 384bits, we need 512bits, or 1.333bits to save a single bit
 */
struct Bitvector {
    std::vector<uint64_t> superblocks;
    std::vector<uint8_t>  blocks;
    std::vector<uint64_t> bits;
    size_t totalSize;

    size_t size() const noexcept {
        return totalSize;
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

        auto maskedBits = (bits[blockId] << (63 - bitId));
        auto bitcount   = std::bitset<64>{maskedBits}.count();

        return superblocks[superblockId]
                + blocks[blockId]
                + bitcount;
    }

    template <typename CB>
    Bitvector(size_t length, CB cb) {
        totalSize = length;
        superblocks.reserve(length/(64*4) + 1);
        blocks.reserve(length/64 + 1);
        bits.reserve(length/64 + 1);

        uint64_t sblock_acc{};
        uint16_t block_acc{};

        for (size_t size{0}; size < length; ++size) {
            if (size % 256 == 0) { // new super block + new block
                superblocks.emplace_back(sblock_acc);
                blocks.emplace_back();
                bits.emplace_back();
                block_acc = 0;
            } else if (size % 64 == 0) { // new block
                blocks.emplace_back(block_acc);
            }

            auto bitId        = size % 64;
            auto blockId      = size / 64;

            if (cb(size)) {
                auto& tbits = bits[blockId];
                tbits = tbits | (1ull << bitId);

                block_acc  += 1;
                sblock_acc += 1;
            }
        }
    }

    Bitvector() {}
    Bitvector(Bitvector const&) = default;
    Bitvector(Bitvector&&) noexcept = default;
    auto operator=(Bitvector const&) -> Bitvector& = default;
    auto operator=(Bitvector&&) noexcept -> Bitvector& = default;

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(superblocks, totalSize);
    }
};

static_assert(BitVector_c<Bitvector>);

}
}
