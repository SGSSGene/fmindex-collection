// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "concepts.h"
#include "utils.h"

#include <bitset>
#include <cassert>
#include <vector>

namespace fmindex_collection::rankvector {

template <size_t TSigma, uint64_t TAlignment, typename block_t>
struct InterleavedEPR {
    static_assert(TSigma > 0, "Alphabet has to have at least 1 letter");
    static constexpr size_t Sigma = TSigma;

    // number of full length bitvectors needed `2^bitct ≥ TSigma`
    static constexpr auto bitct = required_bits(TSigma-1);

    // next full power of 2
    static constexpr auto bvct  = pow(2, bitct);

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


        uint64_t prefix_rank(uint64_t idx, uint8_t symb) const {
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
        void serialize(Archive& ar) {
            ar(blocks, inBlock);
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

    uint64_t rank(uint64_t idx, uint8_t symb) const {
        auto blockId      = idx / letterFit;
        auto superBlockId = idx / block_size;
        auto bitId        = idx % letterFit;
        return blocks[blockId].prefix_rank(bitId, symb)
               - ((symb>0)?(blocks[blockId].prefix_rank(bitId, symb-1)):0)
               + superBlocks[superBlockId][symb];
    }

    uint64_t prefix_rank(uint64_t idx, uint8_t symb) const {
        auto blockId      = idx / letterFit;
        auto superBlockId = idx / block_size;
        auto bitId        = idx % letterFit;
        uint64_t a={};
        for (uint64_t i{0}; i<= symb; ++i) {
            a += superBlocks[superBlockId][i];
        }
        return blocks[blockId].prefix_rank(bitId, symb) + a;
    }


    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
        auto blockId      = idx / letterFit;
        auto superBlockId = idx / block_size;
        auto bitId        = idx % letterFit;
        auto res = std::array<uint64_t, TSigma>{};

        auto pre = std::array<uint64_t, TSigma>{};
        for (uint64_t symb{0}; symb < TSigma; ++symb) {
            pre[symb] = blocks[blockId].prefix_rank(bitId, symb);
        }

        res[0] = pre[0] + superBlocks[superBlockId][0];

        for (uint64_t symb{1}; symb < TSigma; ++symb) {
            res[symb] = pre[symb] - pre[symb-1] + superBlocks[superBlockId][symb];
        }
        return res;
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, TSigma>, std::array<uint64_t, TSigma>> {
        auto blockId      = idx / letterFit;
        auto superBlockId = idx / block_size;
        auto bitId        = idx % letterFit;

        auto rs  = std::array<uint64_t, TSigma>{};
        auto prs = std::array<uint64_t, TSigma>{};

        auto pre = std::array<uint64_t, TSigma>{};
        for (uint64_t symb{0}; symb < TSigma; ++symb) {
            pre[symb] = blocks[blockId].prefix_rank(bitId, symb);
        }

        rs[0] = pre[0];
        for (uint64_t symb{1}; symb < TSigma; ++symb) {
            rs[symb]  = pre[symb]-pre[symb-1];
        }
        for (uint64_t symb{0}; symb < TSigma; ++symb) {
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

template <size_t TSigma> using InterleavedEPR8         = InterleavedEPR<TSigma,  8,  uint8_t>;
template <size_t TSigma> using InterleavedEPR16        = InterleavedEPR<TSigma,  8, uint16_t>;
template <size_t TSigma> using InterleavedEPR32        = InterleavedEPR<TSigma,  8, uint32_t>;
template <size_t TSigma> using InterleavedEPR8Aligned  = InterleavedEPR<TSigma, 64,  uint8_t>;
template <size_t TSigma> using InterleavedEPR16Aligned = InterleavedEPR<TSigma, 64, uint16_t>;
template <size_t TSigma> using InterleavedEPR32Aligned = InterleavedEPR<TSigma, 64, uint32_t>;


static_assert(checkRankVector<InterleavedEPR8>);
static_assert(checkRankVector<InterleavedEPR16>);
static_assert(checkRankVector<InterleavedEPR32>);
static_assert(checkRankVector<InterleavedEPR8Aligned>);
static_assert(checkRankVector<InterleavedEPR16Aligned>);
static_assert(checkRankVector<InterleavedEPR32Aligned>);

}
