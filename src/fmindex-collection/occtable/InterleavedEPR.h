// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "../builtins.h"
#include "concepts.h"
#include "utils.h"

#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <span>
#include <vector>


namespace fmindex_collection {
namespace occtable {
namespace interleavedEPR_impl {

template <uint64_t TSigma, uint64_t TAlignment, typename block_t>
struct Bitvector {

    // number of full length bit vectors needed `2^bitct ≥ TSigma`
    static constexpr auto bitct = required_bits(TSigma-1);
    // next full power of 2
    static constexpr auto bvct  = pow(2, bitct);


    // To select a char at the even/uneven position
    static constexpr uint64_t maskEven = []() {
        uint64_t entries = 64 / bitct;
        auto result = uint64_t{0};
        auto chunkMaskEven = uint64_t{(1ull << bitct)-1ull};
        for (uint64_t i{0}; i < entries; i += 2) {
            result = (result << (bitct*2)) | chunkMaskEven;
        }
        return result;
    }();

    static constexpr uint64_t bitMask = []() {
        uint64_t entries = 64 / bitct;
        auto result = uint64_t{0};
        auto mask = 1<<bitct;
        for (uint64_t i{0}; i < entries; i += 2) {
            result = (result << (bitct*2)) | mask;
        }
        return result;

    }();

    static constexpr std::array<uint64_t, TSigma> rb = []() {
        auto results = std::array<uint64_t, TSigma>{};

        for (uint64_t symb{0}; symb < TSigma; ++symb) {
            uint64_t mask = symb | (1<<bitct);
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
//            __builtin_prefetch((const void*)(&bits), 0, 0);
        }


        uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
            assert(idx < 64 / bitct);

            auto _inblock = inBlock;// & ((1ull<<(idx*bitct)) -1);

            auto te = ((rb[symb] - (_inblock & maskEven)) & bitMask) >> bitct;
            auto to = (rb[symb] - ((_inblock>>bitct) & maskEven)) & bitMask;
            auto epr = (te | to) & ((1ull << (idx*bitct))-1ull);

            auto ct = std::bitset<64>{epr}.count();

            for (uint64_t i{0}; i <= symb; ++i) {
                ct += blocks[i];
            }
            return ct;
        }


        uint64_t symbol(uint64_t idx) const {
            assert(idx < 64 / bitct);

            auto mask = uint64_t{(1ull<<bitct)-1ull};
            uint64_t symb = (inBlock >> (idx*bitct)) & mask;
            return symb;
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar(blocks, inBlock);
        }
    };

    constexpr static uint64_t letterFit = 64 / bitct;
    static constexpr uint64_t block_size = ((1ull<<(sizeof(block_t)*8)) / letterFit)*letterFit;

    std::vector<Block> blocks;
    std::vector<std::array<uint64_t, TSigma>> superBlocks;
    std::array<uint64_t, TSigma+1> C;


    Bitvector(std::span<uint8_t const> _bwt) {
        blocks.reserve(_bwt.size()/block_size+1);

        auto sblock_acc = std::array<uint64_t, TSigma>{}; // accumulator for super blocks
        auto block_acc  = std::array<block_t, TSigma>{};  // accumulator for blocks

        for (uint64_t size{0}; size < _bwt.size();) {
            superBlocks.emplace_back(sblock_acc);
            block_acc = {};

            for (uint64_t blockId{0}; blockId < block_size/letterFit and size < _bwt.size(); ++blockId) {
                blocks.emplace_back();
                blocks.back().blocks = block_acc;

                for (uint64_t bitId{0}; bitId < letterFit and size < _bwt.size(); ++bitId, ++size) {

                    uint64_t symb = _bwt[size];
                    blocks.back().inBlock |= symb << (bitct * bitId);

                    block_acc[symb] += 1;
                    sblock_acc[symb] += 1;
                }
            }
        }
        // Add a new block, so we can access one row more than our bwt length (!TODO this might be more than required)
        superBlocks.emplace_back(sblock_acc);
        blocks.emplace_back();
        blocks.back().blocks = block_acc;

        C[0] = 0;
        for (uint64_t i{0}; i < TSigma; ++i) {
            C[i+1] = sblock_acc[i] + C[i];
        }
    }

    Bitvector(cereal_tag) {}


    uint64_t memoryUsage() const {
        return blocks.size() * sizeof(blocks.back())
            + superBlocks.size() * sizeof(superBlocks.back())
            + sizeof(C);
    }

    void prefetch(uint64_t idx) const {
        auto blockId      = idx / letterFit;
        blocks[blockId].prefetch();
    }

    uint64_t rank(uint64_t idx, uint64_t symb) const {
        auto blockId      = idx / letterFit;
        auto superBlockId = idx / block_size;
        auto bitId        = idx % letterFit;
        return blocks[blockId].prefix_rank(bitId, symb)
               - ((symb>0)?(blocks[blockId].prefix_rank(bitId, symb-1)):0)
               + superBlocks[superBlockId][symb]
               + C[symb];
    }

    uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
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

        res[0] = pre[0] + superBlocks[superBlockId][0] + C[0];

        for (uint64_t symb{1}; symb < TSigma; ++symb) {
            res[symb] = pre[symb] - pre[symb-1] + superBlocks[superBlockId][symb] + C[symb];
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

        for (uint64_t symb{0}; symb < TSigma; ++symb) {
            rs[symb] += C[symb];
        }
        return {rs, prs};
    }

    uint64_t symbol(uint64_t idx) const {
        auto blockId      = idx / letterFit;
        auto bitId        = idx % letterFit;
        return blocks[blockId].symbol(bitId);
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(blocks, superBlocks, C);
    }
};


template <uint64_t TSigma, typename block_t, uint64_t TAlignment>
struct OccTable {
    using TLengthType = uint64_t;
    static constexpr uint64_t Sigma = TSigma;

    Bitvector<Sigma, TAlignment, block_t> bitvector;

    static uint64_t expectedMemoryUsage(uint64_t length) {
        using Block = typename Bitvector<TSigma, TAlignment, block_t>::Block;
        auto blockSize = std::max(alignof(Block), sizeof(Block));

        uint64_t C           = sizeof(uint64_t) * (Sigma+1);
        uint64_t blocks      = blockSize        * (length+1) / 64;
        uint64_t superblocks = sizeof(uint64_t) * (length+1) / (1ull << (sizeof(block_t) * 8));
        return C + blocks + superblocks;
    }

    OccTable(std::span<uint8_t const> _bwt)
        : bitvector{_bwt}
    {}

    OccTable(cereal_tag)
        : bitvector{cereal_tag{}}
    {}

    uint64_t memoryUsage() const {
        return bitvector.memoryUsage() + sizeof(OccTable);
    }

    uint64_t size() const {
        return bitvector.C.back();
    }

    auto prefetch(uint64_t idx) const {
        bitvector.prefetch(idx);
    }

    uint64_t rank(uint64_t idx, uint64_t symb) const {
        return bitvector.rank(idx, symb);
    }

    uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
        return bitvector.prefix_rank(idx, symb);
    }

    uint64_t symbol(uint64_t idx) const {
        return bitvector.symbol(idx);
    }

    auto all_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, Sigma>, std::array<uint64_t, Sigma>> {
        auto [rs, prs] = bitvector.all_ranks_and_prefix_ranks(idx);
        return {rs, prs};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bitvector);
    }
};
}

namespace interleavedEPR8 {
template <uint64_t TSigma>
struct OccTable : interleavedEPR_impl::OccTable<TSigma, uint8_t, 8> {
    using interleavedEPR_impl::OccTable<TSigma, uint8_t, 8>::OccTable;
    static auto name() -> std::string {
        return "Interleaved EPR (8bit)";
    }

    static auto extension() -> std::string {
        return "iepr8";
    }
};
static_assert(checkOccTable<OccTable>);
}
namespace interleavedEPR16 {
template <uint64_t TSigma>
struct OccTable : interleavedEPR_impl::OccTable<TSigma, uint16_t, 8> {
    using interleavedEPR_impl::OccTable<TSigma, uint16_t, 8>::OccTable;

    static auto name() -> std::string {
        return "Interleaved EPR (16bit)";
    }

    static auto extension() -> std::string {
        return "iepr16";
    }
};
static_assert(checkOccTable<OccTable>);
}
namespace interleavedEPR32 {
template <uint64_t TSigma>
struct OccTable : interleavedEPR_impl::OccTable<TSigma, uint32_t, 8> {
    using interleavedEPR_impl::OccTable<TSigma, uint32_t, 8>::OccTable;
    static auto name() -> std::string {
        return "Interleaved EPR (32bit)";
    }

    static auto extension() -> std::string {
        return "iepr32";
    }
};
static_assert(checkOccTable<OccTable>);
}
namespace interleavedEPR8Aligned {
template <uint64_t TSigma>
struct OccTable : interleavedEPR_impl::OccTable<TSigma, uint8_t, 64> {
    using interleavedEPR_impl::OccTable<TSigma, uint8_t, 64>::OccTable;
    static auto name() -> std::string {
        return "Interleaved EPR (8bit, aligned)";
    }

    static auto extension() -> std::string {
        return "iepr8a";
    }
};
static_assert(checkOccTable<OccTable>);
}
namespace interleavedEPR16Aligned {
template <uint64_t TSigma>
struct OccTable : interleavedEPR_impl::OccTable<TSigma, uint16_t, 64> {
    using interleavedEPR_impl::OccTable<TSigma, uint16_t, 64>::OccTable;
    static auto name() -> std::string {
        return "Interleaved EPR (16bit, aligned)";
    }

    static auto extension() -> std::string {
        return "iepr16a";
    }
};
static_assert(checkOccTable<OccTable>);
}
namespace interleavedEPR32Aligned {
template <uint64_t TSigma>
struct OccTable : interleavedEPR_impl::OccTable<TSigma, uint32_t, 64> {
    using interleavedEPR_impl::OccTable<TSigma, uint32_t, 64>::OccTable;
    static auto name() -> std::string {
        return "Interleaved EPR (32bit, aligned)";
    }

    static auto extension() -> std::string {
        return "iepr32a";
    }
};
static_assert(checkOccTable<OccTable>);

}

}
}
