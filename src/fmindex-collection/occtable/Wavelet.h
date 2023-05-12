// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "concepts.h"
#include "utils.h"

#include <array>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <span>
#include <tuple>
#include <vector>

namespace fmindex_collection {
namespace occtable {
namespace wavelet {

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
auto construct_bitvectors(std::span<uint8_t const> bwt) -> std::tuple<std::array<Bitvector, pow(2, required_bits(TSigma))>, std::array<uint64_t, TSigma+1>> {
    constexpr auto bits = required_bits(TSigma);
    constexpr auto bvct = pow(2, bits);
    std::array<Bitvector, bvct> bv;

    auto length = bwt.size();

    auto which_bv = [](uint64_t symb, auto cb) {
        uint64_t id{0};
        uint64_t factor{1};
        uint64_t mask = 1u << (bits-1);
        for (uint64_t b{0}; b < bits; ++b) {
            auto bit = (symb & mask) != 0;
            cb(id, bit);
            id += (bit + 1) * factor;
            factor = factor * 2;
            mask = mask >> 1;
        }
    };

    std::array<uint64_t, TSigma> symb_count{};
    for (uint64_t size{0}; size < length; ++size) {
        auto symb = bwt[size];
        symb_count[symb] += 1;
    }

    std::array<uint64_t, bvct> count{};
    for (uint64_t i{0}; i < TSigma; ++i) {
        which_bv(i, [&](uint64_t id, [[maybe_unused]] uint64_t bit) {
            count[id] += symb_count[i];
        });
    }

    for (uint64_t j{0}; j < bvct; ++j) {
        bv[j].superblocks.reserve(count[j]/384+1);
        bv[j].superblocks.emplace_back();
    }

    std::array<uint64_t, bvct> sblock_acc{0};
    std::array<uint16_t, bvct> block_acc{0};
    std::array<uint64_t, bvct> bv_size{};
    for (auto& c : bv_size) {
        c = 1;
    }

    auto add_bv_bit = [&](uint64_t bv_id, uint64_t value) {
        auto size = bv_size[bv_id];
        if (size % 384 == 0) { // new super block + new block
            bv[bv_id].superblocks.emplace_back();
            bv[bv_id].superblocks.back().superBlockEntry = sblock_acc[bv_id];
            block_acc[bv_id] = 0;
        } else if (size % 64 == 0) { // new block
            bv[bv_id].superblocks.back().setBlock((size % 384) / 64, block_acc[bv_id]);
        }

        auto blockId      = (size >>  6) % 6;
        auto bitId        = size &  63;

        auto& bits = bv[bv_id].superblocks.back().bits[blockId];
        bits = bits | (value << bitId);

        block_acc[bv_id]  += value;
        sblock_acc[bv_id] += value;
        bv_size[bv_id] += 1;
    };

    for (uint64_t size{1}; size <= length; ++size) {
        auto symb = bwt[size-1];
        which_bv(symb, [&](uint64_t id, uint64_t bit) {
            add_bv_bit(id, bit);
        });
    }

    std::array<uint64_t, TSigma+1> C;
    C[0] = 0;
    for (uint64_t i{0}; i < TSigma; ++i) {
        C[i+1] = symb_count[i] + C[i];
    }

    return {std::move(bv), C};
}

template <uint64_t TSigma>
struct OccTable {
    using TLengthType = uint64_t;
    static constexpr uint64_t Sigma = TSigma;
    static constexpr auto bits = required_bits(TSigma);
    static constexpr auto bvct = pow(2, bits);

    std::array<Bitvector, bvct> bitvector;
    std::array<uint64_t, Sigma+1> C;

    static uint64_t expectedMemoryUsage(uint64_t length) {
        uint64_t blockSize   = std::max(alignof(Superblock), sizeof(Superblock));

        uint64_t C           = sizeof(uint64_t) * (Sigma+1);
        uint64_t blocks      = blockSize        * (length+1) / 384 * bits;
        return C + blocks;
    }

    OccTable(std::span<uint8_t const> _bwt) {
        std::tie(bitvector, C) = construct_bitvectors<Sigma>(_bwt);
    }

    OccTable(cereal_tag) {}

    static auto name() -> std::string {
        return "Wavelettrees";
    }

    static auto extension() -> std::string {
        return "wt";
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

    uint64_t rank(uint64_t idx, uint64_t symb) const {
        auto which_bv = [](uint64_t symb, auto cb) {
            uint64_t id{0};
            uint64_t factor{1};
            uint64_t mask = 1u << (bits-1);
            for (uint64_t b{0}; b < bits; ++b) {
                auto bit = (symb & mask) != 0;
                cb(id, bit);
                id += (bit + 1) * factor;
                factor = factor * 2;
                mask = mask >> 1;
            }
        };
        uint64_t a = idx;
        which_bv(symb, [&](uint64_t id, uint64_t value) {
            auto newIdx = bitvector[id].rank(a);
            if (value == 0) {
                a = a - newIdx;
            } else {
                a = newIdx;
            }
        });
        return a + C[symb];
    }

    uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
        auto which_bv = [](uint64_t symb, auto cb) {
            uint64_t id{0};
            uint64_t factor{1};
            uint64_t mask = 1u << (bits-1);
            for (uint64_t b{0}; b < bits; ++b) {
                auto bit = (symb & mask) != 0;
                cb(id, bit);
                id += (bit + 1) * factor;
                factor = factor * 2;
                mask = mask >> 1;
            }
        };
        uint64_t a{};
        uint64_t pos = idx;
        which_bv(symb+1, [&](uint64_t id, uint64_t value) {
            auto newIdx = bitvector[id].rank(pos);
            if (value == 0) {
                pos = pos - newIdx;
            } else {
                a += pos - newIdx;
                pos = newIdx;
            }
        });

        return a;
    }

    uint64_t symbol(uint64_t idx) const {
        idx += 1;
        for (uint64_t i{0}; i < Sigma-1; ++i) {
            if (rank(idx, i) > rank(idx-1, i)) {
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
