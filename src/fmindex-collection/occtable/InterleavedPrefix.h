// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "concepts.h"

#include <algorithm>
#include <array>
#include <bitset>
#include <cstdint>
#include <span>
#include <vector>

namespace fmindex_collection {
namespace occtable {
namespace interleavedPrefix {

template <uint64_t TSigma>
struct Bitvector {
    struct Block {
        std::array<uint32_t, TSigma> blocks{};
        std::array<uint64_t, TSigma> bits{};

        uint64_t rank(uint8_t idx, uint8_t symb) const {
            auto bitset = std::bitset<64>(bits[symb] << (63-idx));
            return blocks[symb] + bitset.count();
        }

        uint64_t prefix_rank(uint8_t idx, uint8_t symb) const {
            uint64_t b = {};
            uint64_t block = 0;
            for (uint64_t i{0}; i <= symb; ++i) {
                b = b | bits[i];
                block += blocks[i];
            }
            auto bitset = std::bitset<64>(b << (63-idx));
            return block + bitset.count();
        }

        uint8_t symbol(uint8_t idx) const {
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

    std::vector<Block> blocks;
    std::vector<std::array<uint64_t, TSigma>> superBlocks;
    std::array<uint64_t, TSigma+1> C;

    uint64_t memoryUsage() const {
        return blocks.size() * sizeof(blocks.back())
            + superBlocks.size() * sizeof(superBlocks.back())
            + sizeof(C);
    }

    uint64_t rank(uint64_t idx, uint8_t symb) const {
        if (symb == 0) {
            return prefix_rank(idx, symb) + C[symb];
        }
        return prefix_rank(idx, symb) - prefix_rank(idx, symb-1) + C[symb];
    }

    uint64_t prefix_rank(uint64_t idx, uint8_t symb) const {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> 32;
        auto bitId        = idx &  63;
        return blocks[blockId].rank(bitId, symb) + superBlocks[superBlockId][symb];
    }

    uint8_t symbol(uint64_t idx) const {
        idx += 1;
        auto blockId      = idx >>  6;
        auto bitId        = idx &  63;
        return blocks[blockId].symbol(bitId);
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(blocks, superBlocks, C);
    }
};

template <uint64_t TSigma>
Bitvector<TSigma> construct_bitvector(std::span<uint8_t const> _bwt) {
    auto const length = _bwt.size();
    Bitvector<TSigma> bitvector;
    bitvector.blocks.reserve(length/64+2);

    bitvector.blocks.emplace_back();
    bitvector.superBlocks.emplace_back();

    auto& bv = bitvector;

    std::array<uint64_t, TSigma> sblock_acc{0};
    std::array<uint32_t, TSigma> block_acc{0};

    for (uint64_t size{1}; size <= length; ++size) {
        if (size % (1ull<<32) == 0) { // new super block + new block
            bv.superBlocks.emplace_back(sblock_acc);
            bv.blocks.emplace_back();
            block_acc = {};
        } else if (size % 64 == 0) { // new block
            bv.blocks.emplace_back();
            bv.blocks.back().blocks = block_acc;
        }
        auto blockId      = size >>  6;
        auto bitId        = size &  63;

        auto start = _bwt[size-1];
        for (uint64_t symb{start}; symb < TSigma; ++symb) {
            auto& bits = bv.blocks[blockId].bits[symb];
            bits = bits | (1ull << bitId);
            block_acc[symb] += 1;
            sblock_acc[symb] += 1;
        }
    }

    bv.C[0] = 0;
    for (uint64_t i{0}; i < TSigma; ++i) {
        bv.C[i+1] = sblock_acc[i];
    }
    return bv;
}

template <uint64_t TSigma>
struct OccTable {
    using TLengthType = uint64_t;
    static constexpr uint64_t Sigma = TSigma;

    Bitvector<Sigma> bitvector;

    static uint64_t expectedMemoryUsage(uint64_t length) {
        using Block = typename Bitvector<TSigma>::Block;
        auto blockSize = std::max(alignof(Block), sizeof(Block));

        uint64_t C           = sizeof(uint64_t) * (Sigma+1);
        uint64_t blocks      = blockSize        * (length+1) / 64;
        uint64_t superblocks = sizeof(uint64_t) * (length+1) / (1ull << 32);
        return C + blocks + superblocks;
    }

    OccTable(std::span<uint8_t const> _bwt) {
        bitvector = construct_bitvector<Sigma>(_bwt);
    }

    OccTable(cereal_tag) {}

    static auto name() -> std::string {
        return "Interleaved Prefixed";
    }

    static auto extension() -> std::string {
        return "ipre";
    }
    uint64_t memoryUsage() const {
        return bitvector.memoryUsage() + sizeof(OccTable);
    }

    uint64_t size() const {
        return bitvector.C.back();
    }

    uint64_t rank(uint64_t idx, uint8_t symb) const {
        return bitvector.rank(idx, symb);
    }

    uint64_t prefix_rank(uint64_t idx, uint8_t symb) const {
        return bitvector.prefix_rank(idx, symb);
    }

    uint64_t symbol(uint64_t idx) const {
        return bitvector.symbol(idx);
    }

    auto all_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, Sigma>, std::array<uint64_t, Sigma>> {
        std::array<uint64_t, Sigma> rs{0};
        std::array<uint64_t, Sigma> prs{0};
        for (uint64_t i{0}; i < Sigma; ++i) {
            rs[i] = rank(idx, i);
            prs[i] = prefix_rank(idx, i);
        }
        return {rs, prs};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bitvector);
    }
};
static_assert(checkOccTable<OccTable>);

}
}
}
