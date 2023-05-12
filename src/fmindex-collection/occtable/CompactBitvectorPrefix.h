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
#include <cassert>
#include <cstdint>
#include <span>
#include <tuple>
#include <vector>

namespace fmindex_collection {
namespace occtable {
namespace compactBitvectorPrefix {

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

    bool value(uint64_t idx) const noexcept {
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

struct Bitvector {
    std::vector<Superblock> superblocks{};

    uint64_t rank(uint64_t idx) const noexcept {
        auto superblockId = idx / 384;
        auto bitId        = idx % 384;
        return superblocks[superblockId].rank(bitId);
    }

    bool value(uint64_t idx) const noexcept {
        idx += 1;
        auto superblockId = idx / 384;
        auto bitId        = idx % 384;
        return superblocks[superblockId].value(bitId);
    }

    uint64_t memoryUsage() const {
        return superblocks.size() * sizeof(superblocks.back());
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(superblocks);
    }
};


template <uint64_t TSigma>
auto construct_bitvectors(std::span<uint8_t const> bwt) -> std::tuple<std::array<Bitvector, TSigma>, std::array<Bitvector, TSigma>, std::array<uint64_t, TSigma+1>> {
    std::array<Bitvector, TSigma> bv1;
    std::array<Bitvector, TSigma> bv2;

    auto length = bwt.size();


    for (uint64_t j{0}; j < TSigma; ++j) {
        bv1[j].superblocks.reserve(length/384+1);
        bv2[j].superblocks.reserve(length/384+1);
        bv1[j].superblocks.emplace_back();
        bv2[j].superblocks.emplace_back();
    }

    std::array<uint64_t, TSigma> sblock_acc{0};
    std::array<uint16_t, TSigma> block_acc{0};

    std::array<uint64_t, TSigma> sblock_acc2{0};
    std::array<uint16_t, TSigma> block_acc2{0};


    for (uint64_t size{1}; size <= length; ++size) {
        if (size % 384 == 0) { // new super block + new block
            for (uint64_t j{0}; j < TSigma; ++j) {
                bv1[j].superblocks.emplace_back();
                bv1[j].superblocks.back().superBlockEntry = sblock_acc[j];
                block_acc[j] = 0;

                bv2[j].superblocks.emplace_back();
                bv2[j].superblocks.back().superBlockEntry = sblock_acc2[j];
                block_acc2[j] = 0;
            }
        } else if (size % 64 == 0) { // new block
            for (uint64_t j{0}; j < TSigma; ++j) {
                bv1[j].superblocks.back().setBlock((size % 384) / 64, block_acc[j]);
                bv2[j].superblocks.back().setBlock((size % 384) / 64, block_acc2[j]);
            }
        }

        auto blockId      = (size >>  6) % 6;
        auto bitId        = size &  63;

        { // update rank
            auto symb = bwt[size-1];
            auto& bits = bv1[symb].superblocks.back().bits[blockId];
            bits = bits | (1ull << bitId);

            block_acc[symb]  += 1;
            sblock_acc[symb] += 1;
        }
        { // update prefix rank
            auto symb = bwt[size-1];
            for (uint64_t j{symb}; j < TSigma; ++j) {
                auto& bits = bv2[j].superblocks.back().bits[blockId];
                bits = bits | (1ull << bitId);

                block_acc2[j]  += 1;
                sblock_acc2[j] += 1;
            }
        }

    }

    std::array<uint64_t, TSigma+1> C;
    C[0] = 0;
    for (uint64_t i{0}; i < TSigma; ++i) {
        C[i+1] = sblock_acc[i] + C[i];
    }
    return {std::move(bv1), std::move(bv2), C};
}

template <uint64_t TSigma>
struct OccTable {
    using TLengthType = uint64_t;
    static constexpr uint64_t Sigma = TSigma;

    std::array<Bitvector, Sigma> bitvector1;
    std::array<Bitvector, Sigma> bitvector2;

    std::array<uint64_t, Sigma+1> C;

    static uint64_t expectedMemoryUsage(uint64_t length) {
        uint64_t blockSize   = std::max(alignof(Superblock), sizeof(Superblock));
        uint64_t C           = sizeof(uint64_t) * (Sigma+1);
        uint64_t blocks      = blockSize        * (length+1) / 384 * Sigma;
        return C + blocks*2;
    }


    OccTable(std::span<uint8_t const> _bwt) {
        std::tie(bitvector1, bitvector2, C) = construct_bitvectors<Sigma>(_bwt);
    }

    OccTable(cereal_tag) {}

    static auto name() -> std::string {
        return "CompactBitvectorPrefix";
    }

    static auto extension() -> std::string {
        return "cbvpre";
    }

    uint64_t memoryUsage() const {
        uint64_t memory{};
        for (auto const& bv : bitvector1) {
            memory += bv.memoryUsage();
        }
        for (auto const& bv : bitvector2) {
            memory += bv.memoryUsage();
        }
        memory += sizeof(OccTable);
        return memory;
    }

    uint64_t size() const {
        return C.back();
    }

    uint64_t rank(uint64_t idx, uint8_t symb) const {
        return bitvector1[symb].rank(idx) + C[symb];
    }

    uint64_t prefix_rank(uint64_t idx, uint8_t symb) const {
        return bitvector2[symb].rank(idx);
    }

    uint64_t symbol(uint64_t idx) const {
        for (uint64_t i{0}; i < Sigma-1; ++i) {
            if (bitvector1[i].value(idx)) {
                return i;
            }
        }
        return Sigma-1;
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
        ar(bitvector1, bitvector2, C);
    }
};
static_assert(checkOccTable<OccTable>);

}
}
}
