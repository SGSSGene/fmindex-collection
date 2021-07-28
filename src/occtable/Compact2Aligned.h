#pragma once

#include "concepts.h"

#include <array>
#include <bitset>
#include <cstdint>
#include <vector>

namespace occtable {
namespace compact2Aligned {

template <size_t TSigma>
struct Bitvector {
    struct alignas(64) Block {
        std::array<uint16_t, TSigma> blocks{};
        std::array<uint64_t, TSigma> bits{};
        uint64_t rank(uint8_t symb, uint8_t idx) const {
            auto bitset = std::bitset<64>(bits[symb] << (63-idx));
            return blocks[symb] + bitset.count();
        }

        uint64_t prefix_rank(uint8_t symb, uint8_t idx) const {
            uint64_t b = {};
            uint64_t block = 0;
            for (size_t i{0}; i <= symb; ++i) {
                b = b | bits[i];
                block += blocks[i];
            }
            auto bitset = std::bitset<64>(b << (63-idx));
            return block + bitset.count();

        }
    };

    std::vector<Block> blocks;
    std::vector<std::array<uint64_t, TSigma>> superBlocks;
    std::array<uint64_t, TSigma+1> C;

    uint64_t rank(uint8_t symb, uint64_t idx) const {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> 16;
        auto bitId        = idx &  63;
        return blocks[blockId].rank(symb, bitId) + superBlocks[superBlockId][symb] + C[symb];
    }

    uint64_t prefix_rank(uint8_t symb, uint64_t idx) const {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> 16;
        auto bitId        = idx &  63;
        uint64_t a={};
        for (size_t i{0}; i<= symb; ++i) {
            a += superBlocks[superBlockId][i];
        }
        return blocks[blockId].prefix_rank(symb, bitId) + a;
    }

};

template <size_t TSigma, typename CB>
Bitvector<TSigma> construct_bitvector(size_t length, CB cb) {
    Bitvector<TSigma> bitvector;
    bitvector.blocks.reserve(length/64+2);

    bitvector.blocks.emplace_back();
    bitvector.superBlocks.emplace_back();

    auto& bv = bitvector;

    std::array<uint64_t, TSigma> sblock_acc{0};
    std::array<uint16_t, TSigma> block_acc{0};

    for (size_t size{1}; size <= length; ++size) {
        if (size % (1ul<<16) == 0) { // new super block + new block
            bv.superBlocks.emplace_back(sblock_acc);
            bv.blocks.emplace_back();
            block_acc = {};
        } else if (size % 64 == 0) { // new block
            bv.blocks.emplace_back();
            bv.blocks.back().blocks = block_acc;
        }
        auto blockId      = size >>  6;
        auto superBlockId = size >> 16;
        auto bitId        = size &  63;

        auto symb = cb(size-1);

        auto& bits = bv.blocks[blockId].bits[symb];
        bits = bits | (1ul << bitId);
        block_acc[symb] += 1;
        sblock_acc[symb] += 1;
    }

    bv.C[0] = 0;
    for (size_t i{0}; i < TSigma; ++i) {
        bv.C[i+1] = sblock_acc[i] + bv.C[i];
    }
    return bv;
};

template <size_t TSigma>
struct OccTable {
    static constexpr size_t Sigma = TSigma;

    Bitvector<Sigma> bitvector;

    static size_t expectedMemoryUsage(size_t length) {
        using Block = typename Bitvector<TSigma>::Block;
        auto blockSize = std::max(alignof(Block), sizeof(Block));

        size_t C           = sizeof(uint64_t) * (Sigma+1);
        size_t blocks      = blockSize        * (length+1) / 64;
        size_t superblocks = sizeof(uint64_t) * (length+1) / (1ul << 16);
        return C + blocks + superblocks;
    }


    OccTable(std::vector<uint8_t> const& _bwt) {
        bitvector = construct_bitvector<Sigma>(_bwt.size(), [&](size_t i) -> uint8_t {
            return _bwt[i];
        });
    }

    uint64_t size() const {
        return bitvector.C.back();
    }

    uint64_t rank(uint8_t symb, uint64_t idx) const {
        return bitvector.rank(symb, idx);
    }

    uint64_t prefix_rank(uint8_t symb, uint64_t idx) const {
        return bitvector.prefix_rank(symb, idx);
    }
};
static_assert(checkOccTable<OccTable>);

}
}
