// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "../builtins.h"
#include "concepts.h"

#include <algorithm>
#include <array>
#include <bitset>
#include <cstdint>
#include <span>
#include <vector>

/**
 * A generic Interleaved OccTable
 *
 * This interleaves blocks and bits. Blocks size is determined by `block_t`.
 *
 * Example dna5+sentinel with block_t=uint32_t
 * superblocks: 64bit
 * interleaved blocks+bits: (32bit * 6) + (64bit * 6)
 *
 *
 */
namespace fmindex_collection {
namespace occtable {
namespace genericInterleaved {

template <uint64_t TSigma, uint64_t TAlignment, typename block_t>
struct Bitvector {
    struct alignas(TAlignment) Block {
        std::array<block_t, TSigma> blocks{};
        std::array<uint64_t, TSigma> bits{};

        void prefetch() const {
            __builtin_prefetch(reinterpret_cast<void const*>(&blocks), 0, 0);
//            __builtin_prefetch((const void*)(&bits), 0, 0);
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

    static constexpr uint64_t block_size = sizeof(block_t) * 8;

    std::vector<Block> blocks;
    std::vector<std::array<uint64_t, TSigma>> superBlocks;
    std::array<uint64_t, TSigma+1> C;

    Bitvector(std::span<uint8_t const> bwt) {
        auto length = bwt.size();
        blocks.reserve(length/64+2);

        blocks.emplace_back();
        superBlocks.emplace_back();

        std::array<uint64_t, TSigma> sblock_acc{0};
        std::array<block_t, TSigma> block_acc{0};

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

            auto symb = bwt[size-1];

            auto& bits = blocks[blockId].bits[symb];
            bits = bits | (1ull << bitId);
            block_acc[symb] += 1;
            sblock_acc[symb] += 1;
        }

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
        auto blockId      = idx >>  6;
        blocks[blockId].prefetch();
    }

    uint64_t rank(uint64_t idx, uint64_t symb) const {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;
        auto bitId        = idx &  63;
        return blocks[blockId].rank(bitId, symb) + superBlocks[superBlockId][symb] + C[symb];
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
            res[symb] = blocks[blockId].rank(bitId, symb) + superBlocks[superBlockId][symb] + C[symb];
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

        for (uint64_t symb{0}; symb < TSigma; ++symb) {
            rs[symb] += C[symb];
        }
        return {rs, prs};
    }

    uint64_t symbol(uint64_t idx) const {
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


template <uint64_t TSigma, uint64_t TAlignment, typename block_t>
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
}
}
