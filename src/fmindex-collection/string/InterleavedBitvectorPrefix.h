// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "concepts.h"

#include <bitset>
#include <vector>

namespace fmc::string {

template <size_t TSigma, uint64_t TAlignment, typename block_t>
struct InterleavedBitvectorPrefix {
    struct alignas(TAlignment) Block {
        std::array<block_t, TSigma> blocks{};
        std::array<uint64_t, TSigma> bits{};

        void prefetch() const {
            __builtin_prefetch(reinterpret_cast<void const*>(&blocks), 0, 0);
        }

        uint64_t rank(uint64_t idx, uint64_t symb) const {
            auto _bits = bits[symb];
            auto block = blocks[symb];
            if (symb > 0) {
                _bits = _bits & ~bits[symb-1];
                block = block - blocks[symb-1];
            }
            auto bitset = std::bitset<64>(_bits << (63-idx));
            return block + bitset.count();
        }

        uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
            auto bitset = std::bitset<64>(bits[symb] << (63-idx));
            return blocks[symb] + bitset.count();
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
        void serialize(this auto&& self, Archive& ar) {
            ar(self.blocks, self.bits);
        }
    };

    static constexpr size_t Sigma = TSigma;

    static constexpr uint64_t block_size = sizeof(block_t) * 8;

    std::vector<Block> blocks;
    std::vector<std::array<uint64_t, TSigma>> superBlocks;
    size_t totalLength{};

    InterleavedBitvectorPrefix() = default;
    InterleavedBitvectorPrefix(std::span<uint8_t const> _symbols) {
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

            auto _symb = _symbols[size-1];

            for (size_t symb{_symb}; symb < Sigma; ++symb) {
                auto& bits = blocks[blockId].bits[symb];
                bits = bits | (1ull << bitId);
                block_acc[symb] += 1;
                sblock_acc[symb] += 1;
            }
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
        auto r = blocks[blockId].rank(bitId, symb) + superBlocks[superBlockId][symb];
        if (symb > 0) {
            r -= superBlocks[superBlockId][symb-1];
        }
        return r;
    }

    uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
        if (symb == 0) return 0;
        symb -= 1;
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;
        auto bitId        = idx &  63;
        return blocks[blockId].prefix_rank(bitId, symb) + superBlocks[superBlockId][symb];
    }

    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
        auto res = std::array<uint64_t, TSigma>{};
        for (uint64_t symb{0}; symb < TSigma; ++symb) {
            res[symb] = rank(idx, symb);
        }
        return res;
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, TSigma>, std::array<uint64_t, TSigma>> {
        auto rs  = std::array<uint64_t, TSigma>{};
        auto prs = std::array<uint64_t, TSigma>{};

        for (uint64_t symb{0}; symb < TSigma; ++symb) {
            rs[symb]  = rank(idx, symb);
            prs[symb] = prefix_rank(idx, symb);
        }
        return {rs, prs};
    }

    template <typename Archive>
    void serialize(this auto&& self, Archive& ar) {
        ar(self.blocks, self.superBlocks, self.totalLength);
    }
};

template <size_t TSigma> using InterleavedBitvectorPrefix8         = InterleavedBitvectorPrefix<TSigma,  8,  uint8_t>;
template <size_t TSigma> using InterleavedBitvectorPrefix16        = InterleavedBitvectorPrefix<TSigma,  8, uint16_t>;
template <size_t TSigma> using InterleavedBitvectorPrefix32        = InterleavedBitvectorPrefix<TSigma,  8, uint32_t>;
template <size_t TSigma> using InterleavedBitvectorPrefix8Aligned  = InterleavedBitvectorPrefix<TSigma, 64,  uint8_t>;
template <size_t TSigma> using InterleavedBitvectorPrefix16Aligned = InterleavedBitvectorPrefix<TSigma, 64, uint16_t>;
template <size_t TSigma> using InterleavedBitvectorPrefix32Aligned = InterleavedBitvectorPrefix<TSigma, 64, uint32_t>;


static_assert(checkString_c<InterleavedBitvectorPrefix8>);
static_assert(checkString_c<InterleavedBitvectorPrefix16>);
static_assert(checkString_c<InterleavedBitvectorPrefix32>);
static_assert(checkString_c<InterleavedBitvectorPrefix8Aligned>);
static_assert(checkString_c<InterleavedBitvectorPrefix16Aligned>);
static_assert(checkString_c<InterleavedBitvectorPrefix32Aligned>);

}
