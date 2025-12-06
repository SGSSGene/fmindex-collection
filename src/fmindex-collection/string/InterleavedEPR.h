// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../builtins.h"
#include "concepts.h"

#include <bit>
#include <bitset>
#include <cassert>
#include <vector>

namespace fmc::string {

template <size_t TSigma, uint64_t TAlignment, typename block_t>
struct InterleavedEPR {
    static_assert(TSigma > 0, "Alphabet has to have at least 1 letter");
    static constexpr size_t Sigma = TSigma;

    // number of full length bitvectors needed `2^bitct ≥ TSigma`
    static constexpr auto bitct = std::bit_width(TSigma-1);

    // next full power of 2
    static constexpr auto bvct  = std::bit_ceil(TSigma);

    // To select a char at the even/uneven position
    static constexpr uint64_t maskEven = []() {
        size_t entries = 64 / bitct;
        auto result = uint64_t{0};
        auto chunkMaskEven = (uint64_t{1} << bitct)-1;
        for (size_t i{0}; i < entries; i += 2) {
            result = (result << (bitct*2)) | chunkMaskEven;
        }
        return result;
    }();

    static constexpr uint64_t bitMask = []() {
        auto entries = size_t{64 / bitct};
        auto result = uint64_t{0};
        auto mask = uint64_t{1}<<bitct;
        for (uint64_t i{0}; i < entries; i += 2) {
            result = (result << (bitct*2)) | mask;
        }
        return result;

    }();

    static constexpr std::array<uint64_t, TSigma> rb = []() {
        auto results = std::array<uint64_t, TSigma>{};

        for (uint64_t symb{0}; symb < TSigma; ++symb) {
            uint64_t mask = symb | (uint64_t{1}<<bitct);
            uint64_t entries = 64 / bitct;
            for (uint64_t i{0}; i < entries; i += 2) {
                results[symb] = (results[symb] << (bitct*2)) | mask;
            }
        }

        return results;
    }();

    struct alignas(TAlignment) Block {
        std::array<block_t, TSigma> blocks{};
        uint64_t inBlock{};

        void prefetch() const {
            __builtin_prefetch(reinterpret_cast<void const*>(&blocks), 0, 0);
        }


        uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
            if (symb == 0) return 0;
            symb -= 1;
            assert(idx < 64 / bitct);

            auto _inblock = inBlock;// & ((1ull<<(idx*bitct)) -1);

            auto te = ((rb[symb] - (_inblock & maskEven)) & bitMask) >> bitct;
            auto to = (rb[symb] - ((_inblock>>bitct) & maskEven)) & bitMask;
            auto epr = (te | to) & ((uint64_t{1} << (idx*bitct))-1);

            auto ct = std::bitset<64>{epr}.count();

            for (uint64_t i{0}; i <= symb; ++i) {
                ct += blocks[i];
            }
            return ct;
        }

        uint8_t symbol(uint64_t idx) const {
            assert(idx < 64 / bitct);

            auto mask = (uint64_t{1}<<bitct)-1;
            uint64_t symb = (inBlock >> (idx*bitct)) & mask;
            return symb;
        }

        template <typename Archive>
        void serialize(this auto&& self, Archive& ar) {
            ar(self.blocks, self.inBlock);
        }
    };

    constexpr static uint64_t letterFit = 64 / bitct;
    static constexpr uint64_t block_size = ((uint64_t{1}<<(sizeof(block_t)*8)) / letterFit)*letterFit;

    std::vector<Block> blocks;
    std::vector<std::array<uint64_t, TSigma>> superBlocks;
    size_t totalLength{};


    InterleavedEPR() = default;
    InterleavedEPR(std::span<uint8_t const> _symbols) {
        totalLength = _symbols.size();
        blocks.reserve(_symbols.size()/block_size+1);

        auto sblock_acc = std::array<uint64_t, TSigma>{}; // accumulator for super blocks
        auto block_acc  = std::array<block_t, TSigma>{};  // accumulator for blocks

        for (uint64_t size{0}; size < _symbols.size();) {
            superBlocks.emplace_back(sblock_acc);
            block_acc = {};

            for (uint64_t blockId{0}; blockId < block_size/letterFit and size < _symbols.size(); ++blockId) {
                blocks.emplace_back();
                blocks.back().blocks = block_acc;

                for (uint64_t bitId{0}; bitId < letterFit and size < _symbols.size(); ++bitId, ++size) {

                    uint64_t symb = _symbols[size];
                    blocks.back().inBlock |= symb << (bitct * bitId);

                    block_acc[symb] += 1;
                    sblock_acc[symb] += 1;
                }
            }
        }
        // Add a new block, so we can access one row more than our symbols length (!TODO this might be more than required)
        superBlocks.emplace_back(sblock_acc);
        blocks.emplace_back();
        blocks.back().blocks = block_acc;
    }

    void prefetch(uint64_t idx) const {
        auto blockId      = idx / letterFit;
        blocks[blockId].prefetch();
    }

    size_t size() const {
        return totalLength;
    }

    uint8_t symbol(uint64_t idx) const {
        auto blockId      = idx / letterFit;
        auto bitId        = idx % letterFit;
        return blocks[blockId].symbol(bitId);
    }

    uint64_t rank(uint64_t idx, uint64_t symb) const {
        auto blockId      = idx / letterFit;
        auto superBlockId = idx / block_size;
        auto bitId        = idx % letterFit;
        return blocks[blockId].prefix_rank(bitId, symb+1)
               - blocks[blockId].prefix_rank(bitId, symb)
               + superBlocks[superBlockId][symb];
    }

    uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
        auto blockId      = idx / letterFit;
        auto superBlockId = idx / block_size;
        auto bitId        = idx % letterFit;
        uint64_t a{};
        for (uint64_t i{0}; i < symb; ++i) {
            a += superBlocks[superBlockId][i];
        }
        return blocks[blockId].prefix_rank(bitId, symb) + a;
    }


    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
        auto blockId      = idx / letterFit;
        auto superBlockId = idx / block_size;
        auto bitId        = idx % letterFit;

        auto prs = std::array<uint64_t, TSigma+1>{};
        for (uint64_t symb{0}; symb < TSigma+1; ++symb) {
            prs[symb] = blocks[blockId].prefix_rank(bitId, symb);
        }

        auto rs = std::array<uint64_t, TSigma>{};
        for (uint64_t symb{0}; symb < TSigma; ++symb) {
            rs[symb] = prs[symb+1] - prs[symb] + superBlocks[superBlockId][symb];
        }
        return rs;
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, TSigma>, std::array<uint64_t, TSigma>> {
        auto blockId      = idx / letterFit;
        auto superBlockId = idx / block_size;
        auto bitId        = idx % letterFit;

        auto prs = std::array<uint64_t, TSigma>{};
        for (uint64_t symb{0}; symb < TSigma; ++symb) {
            prs[symb] = blocks[blockId].prefix_rank(bitId, symb);
        }

        auto rs = std::array<uint64_t, TSigma>{};
        for (uint64_t symb{0}; symb+1 < TSigma; ++symb) {
            rs[symb] = prs[symb+1] - prs[symb] + superBlocks[superBlockId][symb];
        }
        rs[TSigma-1] = blocks[blockId].prefix_rank(bitId, TSigma) - prs[TSigma-1] + superBlocks[superBlockId][TSigma-1];
        return {rs, prs};
    }

    template <typename Archive>
    void serialize(this auto&& self, Archive& ar) {
        ar(self.blocks, self.superBlocks, self.totalLength);
    }
};

template <size_t TSigma> using InterleavedEPR8         = InterleavedEPR<TSigma,  8,  uint8_t>;
template <size_t TSigma> using InterleavedEPR16        = InterleavedEPR<TSigma,  8, uint16_t>;
template <size_t TSigma> using InterleavedEPR32        = InterleavedEPR<TSigma,  8, uint32_t>;
template <size_t TSigma> using InterleavedEPR8Aligned  = InterleavedEPR<TSigma, 64,  uint8_t>;
template <size_t TSigma> using InterleavedEPR16Aligned = InterleavedEPR<TSigma, 64, uint16_t>;
template <size_t TSigma> using InterleavedEPR32Aligned = InterleavedEPR<TSigma, 64, uint32_t>;


static_assert(checkString_c<InterleavedEPR8>);
static_assert(checkString_c<InterleavedEPR16>);
static_assert(checkString_c<InterleavedEPR32>);
static_assert(checkString_c<InterleavedEPR8Aligned>);
static_assert(checkString_c<InterleavedEPR16Aligned>);
static_assert(checkString_c<InterleavedEPR32Aligned>);

}
