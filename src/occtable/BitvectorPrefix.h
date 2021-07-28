#pragma once

#include "concepts.h"

#include <array>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <tuple>
#include <vector>

namespace occtable {
namespace bitvectorPrefix {

struct alignas(64) Superblock {
    uint64_t superBlockEntry{};
    uint64_t blockEntries{};
    std::array<uint64_t, 6> bits{};

    uint64_t rank(size_t idx) const noexcept {
        assert(idx < 384);

        auto blockId = idx >> 6;
        auto block = 0b111111111ul & (blockEntries >> (blockId * 9));
        auto keep = (idx & 63);
        auto maskedBits = bits[blockId] << (63-keep);
        auto ct = std::bitset<64>{maskedBits}.count();

        auto total = superBlockEntry + block + ct;
        return total;
    }

    void setBlock(size_t blockId, size_t value) {
        blockEntries = blockEntries & ~uint64_t{0b111111111ul << blockId*9};
        blockEntries = blockEntries | uint64_t{value << blockId*9};
    }
};

struct Bitvector {
    std::vector<Superblock> superblocks{};

    uint64_t rank(size_t idx) const noexcept {
        auto superblockId = idx / 384;
        auto bitId        = idx % 384;
        return superblocks[superblockId].rank(bitId);
    }
};


template <size_t TSigma, typename CB>
auto construct_bitvectors(size_t length, CB cb) -> std::tuple<std::array<Bitvector, TSigma>, std::array<Bitvector, TSigma>, std::array<uint64_t, TSigma+1>> {
    std::array<Bitvector, TSigma> bv1;
    std::array<Bitvector, TSigma> bv2;


    for (size_t j{0}; j < TSigma; ++j) {
        bv1[j].superblocks.reserve(length/384+1);
        bv2[j].superblocks.reserve(length/384+1);
        bv1[j].superblocks.emplace_back();
        bv2[j].superblocks.emplace_back();
    }

    std::array<uint64_t, TSigma> sblock_acc{0};
    std::array<uint16_t, TSigma> block_acc{0};

    std::array<uint64_t, TSigma> sblock_acc2{0};
    std::array<uint16_t, TSigma> block_acc2{0};


    for (size_t size{1}; size <= length; ++size) {
        if (size % 384 == 0) { // new super block + new block
            for (size_t j{0}; j < TSigma; ++j) {
                bv1[j].superblocks.emplace_back();
                bv1[j].superblocks.back().superBlockEntry = sblock_acc[j];
                block_acc[j] = 0;

                bv2[j].superblocks.emplace_back();
                bv2[j].superblocks.back().superBlockEntry = sblock_acc2[j];
                block_acc2[j] = 0;
            }
        } else if (size % 64 == 0) { // new block
            for (size_t j{0}; j < TSigma; ++j) {
                bv1[j].superblocks.back().setBlock((size % 384) / 64, block_acc[j]);
                bv2[j].superblocks.back().setBlock((size % 384) / 64, block_acc2[j]);
            }
        }

        auto blockId      = (size >>  6) % 6;
        auto superBlockId = size / 384;
        auto bitId        = size &  63;

        { // update rank
            auto symb = cb(size-1);
            auto& bits = bv1[symb].superblocks.back().bits[blockId];
            bits = bits | (1ul << bitId);

            block_acc[symb]  += 1;
            sblock_acc[symb] += 1;
        }
        { // update prefix rank
            auto symb = cb(size-1);
            for (size_t j{symb}; j < TSigma; ++j) {
                auto& bits = bv2[j].superblocks.back().bits[blockId];
                bits = bits | (1ul << bitId);

                block_acc2[j]  += 1;
                sblock_acc2[j] += 1;
            }
        }

    }

    std::array<uint64_t, TSigma+1> C;
    C[0] = 0;
    for (size_t i{0}; i < TSigma; ++i) {
        C[i+1] = sblock_acc[i] + C[i];
    }
    return {std::move(bv1), std::move(bv2), C};
};

template <size_t TSigma>
struct OccTable {
    static constexpr size_t Sigma = TSigma;

    std::array<Bitvector, Sigma> bitvector1;
    std::array<Bitvector, Sigma> bitvector2;

    std::array<uint64_t, Sigma+1> C;

    static size_t expectedMemoryUsage(size_t length) {
        size_t blockSize   = std::max(alignof(Superblock), sizeof(Superblock));
        size_t C           = sizeof(uint64_t) * (Sigma+1);
        size_t blocks      = blockSize        * (length+1) / 384 * Sigma;
        return C + blocks*2;
    }


    OccTable(std::vector<uint8_t> const& _bwt) {
        std::tie(bitvector1, bitvector2, C) = construct_bitvectors<Sigma>(_bwt.size(), [&](size_t i) -> uint8_t {
            return _bwt[i];
        });
    }

    uint64_t size() const {
        return C.back();
    }

    uint64_t rank(uint8_t symb, uint64_t idx) const {
        return bitvector1[symb].rank(idx) + C[symb];
    }

    uint64_t prefix_rank(uint8_t symb, uint64_t idx) const {
        return bitvector2[symb].rank(idx);
    }
};
static_assert(checkOccTable<OccTable>);

}
}
