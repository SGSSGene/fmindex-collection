#pragma once

#include "concepts.h"

#include <array>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <tuple>
#include <vector>

namespace occtable {
namespace wavelet {

constexpr inline size_t bits_count(size_t y) {
    if (y == 0) return 1;
    size_t i{0};
    while (y != 0) {
        y = y >> 1;
        ++i;
    }
    return i;
}
constexpr inline size_t pow(size_t b, size_t y) {
    if (y == 0) return 1;
    return pow(b, (y-1)) * b;
}


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

    bool value(size_t idx) const noexcept {
        assert(idx < 384);

        auto blockId = idx >> 6;
        auto bitId = idx & 63;
        return bits[blockId] & (1ul << bitId);
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

    bool value(size_t idx) const noexcept {
        idx += 1;
        auto superblockId = idx / 384;
        auto bitId        = idx % 384;
        return superblocks[superblockId].value(bitId);
    }

    size_t memoryUsage() const {
        return superblocks.size() * sizeof(superblocks.back());
    }
};


template <size_t TSigma, typename CB>
auto construct_bitvectors(size_t length, CB cb) -> std::tuple<std::array<Bitvector, pow(2, bits_count(TSigma))>, std::array<uint64_t, TSigma+1>> {
    constexpr auto bits = bits_count(TSigma);
    constexpr auto bvct = pow(2, bits);
    std::array<Bitvector, bvct> bv;

    auto which_bv = [](uint8_t symb, auto cb) {
        size_t id{0};
        size_t factor{1};
        uint8_t mask = 1u << (bits-1);
        for (size_t b{0}; b < bits; ++b) {
            auto bit = (symb & mask) != 0;
            cb(id, bit);
            id += (bit + 1) * factor;
            factor = factor * 2;
            mask = mask >> 1;
        }
    };

    std::array<size_t, TSigma> symb_count{};
    for (size_t size{0}; size < length; ++size) {
        auto symb = cb(size);
        symb_count[symb] += 1;
    }

    std::array<size_t, bvct> count{};
    for (size_t i{0}; i < TSigma; ++i) {
        which_bv(i, [&](size_t id, size_t bit) {
            count[id] += symb_count[i];
        });
    }


    for (size_t j{0}; j < bvct; ++j) {
        bv[j].superblocks.reserve(count[j]/384+1);
        bv[j].superblocks.emplace_back();
    }

    std::array<uint64_t, bvct> sblock_acc{0};
    std::array<uint16_t, bvct> block_acc{0};
    std::array<uint64_t, bvct> bv_size{};
    for (auto& c : bv_size) {
        c = 1;
    }

    auto add_bv_bit = [&](size_t bv_id, size_t value) {
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

    for (size_t size{1}; size <= length; ++size) {
        auto symb = cb(size-1);
        which_bv(symb, [&](size_t id, size_t bit) {
            add_bv_bit(id, bit);
        });
    }

    std::array<uint64_t, TSigma+1> C;
    C[0] = 0;
    for (size_t i{0}; i < TSigma; ++i) {
        C[i+1] = symb_count[i] + C[i];
    }

    return {std::move(bv), C};
};

template <size_t TSigma>
struct OccTable {
    static constexpr size_t Sigma = TSigma;
    static constexpr auto bits = bits_count(TSigma);
    static constexpr auto bvct = pow(2, bits);


    std::array<Bitvector, bvct> bitvector;
    std::array<uint64_t, Sigma+1> C;

    static size_t expectedMemoryUsage(size_t length) {
        size_t blockSize   = std::max(alignof(Superblock), sizeof(Superblock));

        size_t C           = sizeof(uint64_t) * (Sigma+1);
        size_t blocks      = blockSize        * (length+1) / 384 * bits;
        return C + blocks;
    }

    OccTable(std::vector<uint8_t> const& _bwt) {
        std::tie(bitvector, C) = construct_bitvectors<Sigma>(_bwt.size(), [&](size_t i) -> uint8_t {
            return _bwt[i];
        });
    }
    size_t memoryUsage() const {
        size_t memory{};
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
        auto which_bv = [](uint8_t symb, auto cb) {
            size_t id{0};
            size_t factor{1};
            uint8_t mask = 1u << (bits-1);
            for (size_t b{0}; b < bits; ++b) {
                auto bit = (symb & mask) != 0;
                cb(id, bit);
                id += (bit + 1) * factor;
                factor = factor * 2;
                mask = mask >> 1;
            }
        };
        uint64_t a = idx;
        which_bv(symb, [&](size_t id, size_t value) {
            auto newIdx = bitvector[id].rank(a);
            if (value == 0) {
                a = a - newIdx;
            } else {
                a = newIdx;
            }
        });
        return a + C[symb];
    }

    uint64_t prefix_rank(uint64_t idx, uint8_t symb) const {
        auto which_bv = [](uint8_t symb, auto cb) {
            size_t id{0};
            size_t factor{1};
            uint8_t mask = 1u << (bits-1);
            for (size_t b{0}; b < bits; ++b) {
                auto bit = (symb & mask) != 0;
                cb(id, bit);
                id += (bit + 1) * factor;
                factor = factor * 2;
                mask = mask >> 1;
            }
        };
        uint64_t a{};
        uint64_t pos = idx;
        which_bv(symb+1, [&](size_t id, size_t value) {
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

    uint8_t symbol(uint64_t idx) const {
        idx += 1;
        for (size_t i{0}; i < Sigma-1; ++i) {
            if (rank(idx, i) > rank(idx-1, i)) {
                return i;
            }
        }
        return Sigma-1;
    }

    auto all_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, Sigma>, std::array<uint64_t, Sigma>> {
        std::array<uint64_t, Sigma> rs{0};
        std::array<uint64_t, Sigma> prs{0};
        for (size_t i{0}; i < Sigma; ++i) {
            rs[i] = rank(idx, i);
            prs[i] = prefix_rank(idx, i);
        }
        return {rs, prs};
    }

};

static_assert(checkOccTable<OccTable>);

}
}
