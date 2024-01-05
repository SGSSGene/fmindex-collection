// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../cereal_tag.h"
#include "concepts.h"

#include <bitset>
#include <cassert>
#include <vector>

namespace fmindex_collection {
namespace bitvector {

struct CompactBitvector {
    struct alignas(64) Superblock {
        uint64_t superBlockEntry{};
        uint64_t blockEntries{};
        std::array<uint64_t, 6> bits{};

        uint64_t rank(uint64_t idx) const noexcept {
            assert(idx < 384);

            auto blockId = idx >> 6;
            auto block = 0b111111111ull & (blockEntries >> (blockId * 9));
            auto keep = (idx & 63);
            auto maskedBits = bits[blockId] << (63-keep);
            auto ct = std::bitset<64>{maskedBits}.count();

            auto total = superBlockEntry + block + ct;
            return total;
        }

        bool symbol(uint64_t idx) const noexcept {
            assert(idx < 384);

            auto blockId = idx >> 6;
            auto bitId = idx & 63;
            return bits[blockId] & (1ull << bitId);
        }

        void setBlock(uint64_t blockId, uint64_t value) {
            blockEntries = blockEntries & ~uint64_t{0b111111111ull << blockId*9};
            blockEntries = blockEntries | uint64_t{value << blockId*9};
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar(superBlockEntry, blockEntries, bits);
        }
    };

    using TLengthType = uint64_t;
    static constexpr uint64_t Sigma = 2;

    std::vector<Superblock> superblocks{};
    size_t                  totalLength;

    template <typename CB>
    CompactBitvector(size_t length, CB cb) {
        totalLength = length;

        superblocks.reserve(length/384+1);
        superblocks.emplace_back();

        uint64_t sblock_acc{};
        uint64_t block_acc{};

        for (uint64_t size{1}; size <= length; ++size) {
            if (size % 384 == 0) { // new super block + new block
                superblocks.emplace_back();
                superblocks.back().superBlockEntry = sblock_acc;
                block_acc = 0;
            } else if (size % 64 == 0) { // new block
                superblocks.back().setBlock((size % 384) / 64, block_acc);
            }

            auto blockId      = (size >>  6) % 6;
            auto bitId        = size &  63;

            auto sym = cb(size-1);

            if (sym) {
                auto& bits = superblocks.back().bits[blockId];
                bits = bits | (1ull << bitId);
            }

            block_acc  += 1;
            sblock_acc += 1;
        }
    }

    CompactBitvector(cereal_tag) {}

    CompactBitvector() {}
    CompactBitvector(CompactBitvector const&) = default;
    CompactBitvector(CompactBitvector&&) noexcept = default;
    auto operator=(CompactBitvector const&) -> CompactBitvector& = default;
    auto operator=(CompactBitvector&&) noexcept -> CompactBitvector& = default;


    size_t size() const noexcept {
        return totalLength;
    }

    bool symbol(uint64_t idx) const noexcept {
        idx += 1;
        auto superblockId = idx / 384;
        auto bitId        = idx % 384;
        return superblocks[superblockId].symbol(bitId);
    }

    uint64_t rank(uint64_t idx) const noexcept {
        auto superblockId = idx / 384;
        auto bitId        = idx % 384;
        auto v = superblocks[superblockId].rank(bitId);
        return v;
    }



    template <typename Archive>
    void serialize(Archive& ar) {
        ar(totalLength, superblocks);
    }
};

static_assert(BitVector_c<CompactBitvector>);

}
}
