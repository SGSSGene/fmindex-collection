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
#include <tuple>
#include <span>
#include <vector>

namespace fmindex_collection {
namespace occtable {
namespace bitvector {

struct Bitvector {
    std::vector<uint64_t> superBlockEntry;
    std::vector<uint8_t>  blockEntries;
    std::vector<uint64_t> bits;

    size_t memoryUsage() const {
        return superBlockEntry.size() * 64
        + blockEntries.size() * 8
        + bits.size() * 64;
    }

    uint64_t rank(uint64_t idx) const noexcept {
        auto superblockId = idx / 256;
        auto blockId      = idx / 64;
        auto bitId        = idx & 63;

        auto block = bits[blockId] << (63-bitId);

        auto c = superBlockEntry[superblockId] + blockEntries[blockId]
            + std::bitset<64>{block}.count();
        return c;
    }

    bool value(uint64_t idx) const noexcept {
        idx += 1;

        auto blockId      = idx / 64;
        auto bitId        = idx & 63;
        return bits[blockId] & (1ull << bitId);
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(superBlockEntry, blockEntries, bits);
    }
};

template <uint64_t TSigma>
auto construct_bitvectors(std::span<uint8_t const> bwt) -> std::tuple<std::array<Bitvector, TSigma>, std::array<uint64_t, TSigma+1>> {
    std::array<Bitvector, TSigma> bv;

    auto length = bwt.size();

    for (uint64_t j{0}; j < TSigma; ++j) {
        bv[j].superBlockEntry.reserve(length/256+1);
        bv[j].blockEntries.reserve(length/64+1);
        bv[j].bits.reserve(length/64+1);

        bv[j].superBlockEntry.emplace_back();
        bv[j].blockEntries.emplace_back();
        bv[j].bits.emplace_back();
    }

    std::array<uint64_t, TSigma> sblock_acc{0};
    std::array<uint8_t, TSigma> block_acc{0};

    for (uint64_t size{1}; size <= length; ++size) {
        if (size % 256 == 0) { // new super block + new block
            for (uint64_t j{0}; j < TSigma; ++j) {
                bv[j].superBlockEntry.emplace_back(sblock_acc[j]);
                bv[j].blockEntries.emplace_back();
                bv[j].bits.emplace_back();
                block_acc[j] = 0;
            }
        } else if (size % 64 == 0) { // new block
            for (uint64_t j{0}; j < TSigma; ++j) {
                bv[j].blockEntries.emplace_back(block_acc[j]);
                bv[j].bits.emplace_back();
            }
        }

        auto blockId      = size / 64;;
        auto bitId        = size & 63;

        auto symb = bwt[size-1];
        auto& bits = bv[symb].bits[blockId];
        bits = bits | (1ull << bitId);

        block_acc[symb]  += 1;
        sblock_acc[symb] += 1;
    }

    std::array<uint64_t, TSigma+1> C;
    C[0] = 0;
    for (uint64_t i{0}; i < TSigma; ++i) {
        C[i+1] = sblock_acc[i] + C[i];
    }
    return {std::move(bv), C};
}

template <uint64_t TSigma>
struct OccTable {
    using TLengthType = uint64_t;
    static constexpr uint64_t Sigma = TSigma;

    std::array<Bitvector, Sigma> bitvector;
    std::array<uint64_t, Sigma+1> C;

    static uint64_t expectedMemoryUsage(uint64_t length) {
        uint64_t blockSize   = 256+8*4+64;

        uint64_t C           = sizeof(uint64_t) * (Sigma+1);
        uint64_t blocks      = blockSize        * (length+1) / 256 * Sigma;
        return C + blocks;
    }

    OccTable(std::span<uint8_t const> _bwt) {
        std::tie(bitvector, C) = construct_bitvectors<Sigma>(_bwt);
    }

    OccTable(cereal_tag) {}

    static auto name() -> std::string {
        return "Bitvector";
    }

    static auto extension() -> std::string {
        return "bv";
    }

    uint64_t memoryUsage() const {
        uint64_t memory{};
        for (auto const& bv : bitvector) {
            memory += bv.memoryUsage();
        }
        memory += sizeof(OccTable);
        return memory;
    }

    uint64_t size() const {
        return C.back();
    }

    uint64_t rank(uint64_t idx, uint8_t symb) const {
        return bitvector[symb].rank(idx) + C[symb];
    }

    uint64_t prefix_rank(uint64_t idx, uint8_t symb) const {
        uint64_t a{};
        for (uint64_t i{0}; i <= symb; ++i) {
            a += bitvector[i].rank(idx);
        }
        return a;
    }

    uint64_t symbol(uint64_t idx) const {
        for (uint64_t i{0}; i < Sigma-1; ++i) {
            if (bitvector[i].value(idx)) {
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
        ar(bitvector, C);
    }
};
static_assert(checkOccTable<OccTable>);

}
}
}
