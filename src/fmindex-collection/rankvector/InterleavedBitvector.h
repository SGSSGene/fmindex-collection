// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "concepts.h"

#include <bitset>
#include <vector>

namespace fmindex_collection::rankvector {

template <size_t TSigma, uint64_t TAlignment, typename block_t>
struct InterleavedBitvector {
    struct alignas(TAlignment) Block {
        std::array<block_t, TSigma> blocks{};
        std::array<uint64_t, TSigma> bits{};

        void prefetch() const {
            __builtin_prefetch(reinterpret_cast<void const*>(&blocks), 0, 0);
        }

        uint64_t rank(uint64_t idx, uint64_t symb) const {
            auto bitset = std::bitset<64>(bits[symb] << (63-idx));
            return blocks[symb] + bitset.count();
        }

        uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
            uint64_t b = {};
            uint64_t block = 0;
            for (uint64_t i{0}; i <= symb; ++i) {
                b = b | bits[i];
                block += blocks[i];
            }
            auto bitset = std::bitset<64>(b << (63-idx));
            return block + bitset.count();
        }

        uint64_t symbol(uint64_t idx) const {
            auto bit = (1ull << idx);
            for (uint64_t symb{0}; symb < TSigma-1; ++symb) {
                if (bits[symb] & bit) {
                    return symb;
                }
            }
            return TSigma-1;
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar(blocks, bits);
        }
    };

    static constexpr size_t Sigma = TSigma;

    static constexpr uint64_t block_size = sizeof(block_t) * 8;

    std::vector<Block> blocks;
    std::vector<std::array<uint64_t, TSigma>> superBlocks;
    size_t totalLength;

    InterleavedBitvector() = default;
    InterleavedBitvector(std::span<uint8_t const> _symbols) {
        totalLength = _symbols.size();
        auto length = _symbols.size();
        blocks.reserve(length/64+2);

        blocks.emplace_back();
        superBlocks.emplace_back();

        auto sblock_acc = std::array<uint64_t, TSigma>{0};
        auto block_acc = std::array<block_t, TSigma>{0};

        for (uint64_t size{1}; size <= length; ++size) {
            if (size % (1ull<<block_size) == 0) { // new super block + new block
                superBlocks.emplace_back(sblock_acc);
                blocks.emplace_back();
                block_acc = {};
            } else if (size % 64 == 0) { // new block
                blocks.emplace_back();
                blocks.back().blocks = block_acc;
            }
            auto blockId      = size >>  6;
            auto bitId        = size &  63;

            auto symb = _symbols[size-1];

            auto& bits = blocks[blockId].bits[symb];
            bits = bits | (1ull << bitId);
            block_acc[symb] += 1;
            sblock_acc[symb] += 1;
        }
    }

    size_t size() const noexcept {
        return totalLength;
    }

    void prefetch(uint64_t idx) const {
        auto blockId      = idx >>  6;
        blocks[blockId].prefetch();
    }

    uint8_t symbol(uint64_t idx) const {
        idx += 1;
        auto blockId      = idx >>  6;
        auto bitId        = idx &  63;
        return blocks[blockId].symbol(bitId);
    }

    uint64_t rank(uint64_t idx, uint64_t symb) const {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;
        auto bitId        = idx &  63;
        return blocks[blockId].rank(bitId, symb) + superBlocks[superBlockId][symb];
    }

    uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;
        auto bitId        = idx &  63;
        uint64_t a={};
        for (uint64_t i{0}; i<= symb; ++i) {
            a += superBlocks[superBlockId][i];
        }
        return blocks[blockId].prefix_rank(bitId, symb) + a;
    }

    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;
        auto bitId        = idx &  63;
        auto res = std::array<uint64_t, TSigma>{};
        for (uint64_t symb{0}; symb < TSigma; ++symb) {
            res[symb] = blocks[blockId].rank(bitId, symb) + superBlocks[superBlockId][symb];
        }
        return res;
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, TSigma>, std::array<uint64_t, TSigma>> {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;
        auto bitId        = idx &  63;

        auto rs  = std::array<uint64_t, TSigma>{};
        auto prs = std::array<uint64_t, TSigma>{};

        for (uint64_t symb{0}; symb < TSigma; ++symb) {
            rs[symb]  = blocks[blockId].rank(bitId, symb);
        }
        rs[0] += superBlocks[superBlockId][0];
        for (uint64_t symb{1}; symb < TSigma; ++symb) {
            rs[symb] += superBlocks[superBlockId][symb];
        }
        prs[0] = rs[0];
        for (uint64_t symb{1}; symb < TSigma; ++symb) {
            prs[symb]= prs[symb-1] + rs[symb];
        }
        return {rs, prs};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(blocks, superBlocks, totalLength);
    }
};

template <size_t TSigma> using InterleavedBitvector8         = InterleavedBitvector<TSigma,  8,  uint8_t>;
template <size_t TSigma> using InterleavedBitvector16        = InterleavedBitvector<TSigma,  8, uint16_t>;
template <size_t TSigma> using InterleavedBitvector32        = InterleavedBitvector<TSigma,  8, uint32_t>;
template <size_t TSigma> using InterleavedBitvector8Aligned  = InterleavedBitvector<TSigma, 64,  uint8_t>;
template <size_t TSigma> using InterleavedBitvector16Aligned = InterleavedBitvector<TSigma, 64, uint16_t>;
template <size_t TSigma> using InterleavedBitvector32Aligned = InterleavedBitvector<TSigma, 64, uint32_t>;


static_assert(checkRankVector<InterleavedBitvector8>);
static_assert(checkRankVector<InterleavedBitvector16>);
static_assert(checkRankVector<InterleavedBitvector32>);
static_assert(checkRankVector<InterleavedBitvector8Aligned>);
static_assert(checkRankVector<InterleavedBitvector16Aligned>);
static_assert(checkRankVector<InterleavedBitvector32Aligned>);

}
