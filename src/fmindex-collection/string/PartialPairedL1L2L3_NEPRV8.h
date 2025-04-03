// SPDX-FileCopyrightText: 2024 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../bitset_popcount.h"
#include "../builtins.h"
#include "concepts.h"
#include "utils.h"
#include "EPRV3.h"
#include "NEPRV8.h"

#include <bit>
#include <vector>

#if __has_include(<cereal/types/bitset.hpp>)
#include <cereal/types/bitset.hpp>
#endif

#include <iostream>

#include "../ternarylogic.h"

namespace fmindex_collection::string {

/*
 * PartialPairedL1L2L3 EPRv2 Implementation
 */
template <size_t TSigma, size_t popcount_width, typename blockL0_t = uint16_t, typename blockL1_t = uint32_t>
struct PartialPairedL1L2L3_NEPRV8 {

    static constexpr size_t Sigma = TSigma;

    // number of full length bit vectors needed `2^bitct > TSigma`
    static constexpr auto bitct = std::bit_width(TSigma-1);
    // next full power of 2
    static constexpr auto bvct  = (1ull << bitct);

    static constexpr auto popcount_width_bits = std::bit_width(popcount_width-1);

    struct InBits {
        std::array<std::bitset<popcount_width>, bitct> bits;

        uint8_t symbol(uint64_t idx) const {
            uint8_t symb{};
            for (uint64_t i{bitct}; i > 0; --i) {
                auto b = bits[i-1].test(idx);
                symb = (symb<<1) | b;
            }
            return symb;
        }

        uint64_t rank(uint64_t idx, uint8_t symb) const {
            assert(idx < popcount_width*2);
            auto v = mark_exact_large(symb, bits);
            return signed_rshift_and_count(v, idx);
        }

        uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
            assert(idx < popcount_width*2);

            auto v = mark_exact_or_less_large(symb, bits);
            return signed_rshift_and_count(v, idx);
        }

        auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
            assert(idx < popcount_width*2);

            auto bits2 = bits;
            for (auto& b : bits2) {
                b = b & signed_rightshift_masks<popcount_width>[idx];
            }

            auto vs = neprv8_detail::rank_all<TSigma>(bits2);
            auto v = std::array<uint64_t, TSigma>{};
            static_assert(v.size() <= vs.size());
            for (size_t i{0}; i < v.size(); ++i) {
                v[i] = vs[i].count();
            }
            return v;
        }

        template <typename Archive>
        void load(Archive& ar) {
            for (auto& v : bits) {
                (void)v;
                loadBV(v, ar);
            }
        }
        template <typename Archive>
        void save(Archive& ar) const {
            (void)ar;
            for (auto const& v : bits) {
                (void)v;
                saveBV(v, ar);
            }
        }

    };

    static constexpr uint64_t level0_size = sizeof(blockL0_t) * 8;
    static constexpr uint64_t level1_size = sizeof(blockL1_t) * 8;

    using BlockL0 = std::array<blockL0_t, TSigma>;
    using BlockL1 = std::array<blockL1_t, TSigma>;

    std::vector<InBits> bits;
    std::vector<BlockL0> level0;
    std::vector<BlockL1> level1;
    std::vector<std::array<uint64_t, TSigma>> superBlocks;

    std::vector<BlockL0> level0_prefix;
    std::vector<BlockL1> level1_prefix;
    std::vector<std::array<uint64_t, TSigma>> superBlocks_prefix;


    size_t totalLength;

    PartialPairedL1L2L3_NEPRV8() = default;
    PartialPairedL1L2L3_NEPRV8(std::span<uint8_t const> _symbols) {
        totalLength = _symbols.size();
        auto const length = _symbols.size();
        level1.reserve(length/(1ull<<level1_size)+2);
        level0.reserve(length/popcount_width+2);
        bits.reserve(length/popcount_width+2);

//        auto blockL0_acc = std::array<blockL0_t, TSigma>{0};
        auto blockL1_acc = std::array<blockL1_t, TSigma>{0};
        auto sblock_acc  = std::array<uint64_t, TSigma>{0};

        for (uint64_t size{0}; size < length+1; ++size) {
            if (size % (1ull<<level1_size) == 0) { // new l3 block
                superBlocks.emplace_back(sblock_acc);
                level1.emplace_back();
//                level0.emplace_back();
                if (level0.empty()) {
                    level0.emplace_back();
                }
                level0.back() = {};
                level0.back()[0] = popcount_width;
                bits.emplace_back();
//                blockL0_acc = {};
                blockL1_acc = {};
            } else if (size % (1ull<<level0_size) == 0) { // new l1 block
                level1.emplace_back(blockL1_acc);
//                level0.emplace_back();
                level0.back() = {};
                level0.back()[0] = popcount_width;
                bits.emplace_back();

//                blockL0_acc = {};
            } else if (size % popcount_width == 0) { // new l0 block
                if ((size / popcount_width) % 2 == 1) {
                    level0.emplace_back(level0.back());
                    level0.back()[0] += popcount_width*2;
                }
                bits.emplace_back();
            }
            // Abort, we only wanted to add a new block to the end, if required
            if (size == length) continue;

            auto level0Id     = size >>  (popcount_width_bits);
            auto bitId        = size &  (popcount_width-1);

            uint64_t symb = _symbols[size];

            for (uint64_t i{}; i < bitct; ++i) {
                auto b = ((symb>>i)&1);
                bits[level0Id].bits[i].set(bitId, b);
            }
//            blockL0_acc[symb] += 1;
            level0.back()[symb] += 1;
            blockL1_acc[symb] += 1;
            sblock_acc[symb] += 1;
            level0.back()[0] -= 1;

        }

        // Compute prefix accumulators
        level0_prefix = level0;
        level1_prefix = level1;
        superBlocks_prefix = superBlocks;

        for (auto& vs : level0_prefix) {
            for (size_t i{1}; i < vs.size(); ++i) {
                vs[i] += vs[i-1];
            }
        }
        for (auto& vs : level1_prefix) {
            for (size_t i{1}; i < vs.size(); ++i) {
                vs[i] += vs[i-1];
            }
        }
        for (auto& vs : superBlocks_prefix) {
            for (size_t i{1}; i < vs.size(); ++i) {
                vs[i] += vs[i-1];
            }
        }
    }

    void prefetch(uint64_t idx) const {
        (void)idx;
/*        auto level0Id     = idx >>  6;
        auto level1Id     = idx >> level0_size;
        auto superBlockId = idx >> level1_size;

        __builtin_prefetch(reinterpret_cast<void const*>(&level1[level1Id]), 0, 0);
        __builtin_prefetch(reinterpret_cast<void const*>(&level0[level0Id]), 0, 0);
        __builtin_prefetch(reinterpret_cast<void const*>(&bits[level1Id]), 0, 0);
        __builtin_prefetch(reinterpret_cast<void const*>(&superBlocks[superBlockId]), 0, 0);*/
    }

    size_t size() const {
        return totalLength;
    }

    uint8_t symbol(uint64_t idx) const {
        auto level0Id     = idx >>  popcount_width_bits;
        auto bitId        = idx &  (popcount_width-1);
        return bits[level0Id].symbol(bitId);
    }

    uint64_t rank(uint64_t idx, uint8_t symb) const {
        prefetch(idx);

        auto level0Id     = idx >> popcount_width_bits;
        auto level1Id     = idx >> level0_size;
        auto superBlockId = idx >> level1_size;
        auto bitId        = idx &  (popcount_width*2 - 1);
        auto sign         = (level0Id % 2)*2-1;
        return  size_t{}
                + bits[level0Id].rank(bitId, symb)*sign
                + level0[level0Id/2][symb]
                + level1[level1Id][symb]
                + superBlocks[superBlockId][symb];
    }

    uint64_t prefix_rank(uint64_t idx, uint8_t symb) const {
        prefetch(idx);

        auto level0Id     = idx >>  popcount_width_bits;
        auto level1Id     = idx >> level0_size;
        auto superBlockId = idx >> level1_size;
        auto bitId        = idx &  (popcount_width*2 - 1);
        auto sign         = (level0Id % 2)*2-1;
        return  size_t{}
                + bits[level0Id].prefix_rank(bitId, symb)*sign
                + level0_prefix[level0Id/2][symb]
                + level1_prefix[level1Id][symb]
                + superBlocks_prefix[superBlockId][symb];
    }


    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
        prefetch(idx);

        auto level0Id     = idx >>  popcount_width_bits;
        auto level1Id     = idx >> level0_size;
        auto superBlockId = idx >> level1_size;
        auto bitId        = idx &  (popcount_width*2 - 1);
        auto sign         = (level0Id % 2)*2-1;
        auto res = bits[level0Id].all_ranks(bitId);
        for (uint64_t symb{0}; symb < TSigma; ++symb) {
            res[symb] =   level0[level0Id/2][symb]
                        + level1[level1Id][symb]
                        + superBlocks[superBlockId][symb]
                        + res[symb]*sign;
        }

        return res;
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, TSigma>, std::array<uint64_t, TSigma>> {
        prefetch(idx);

        auto level0Id     = idx >>  popcount_width_bits;
        auto level1Id     = idx >> level0_size;
        auto superBlockId = idx >> level1_size;

        auto prs = std::array<uint64_t, TSigma>{};
        auto bitId        = idx &  (popcount_width*2 - 1);
        auto sign         = (level0Id % 2)*2-1;
        auto res = bits[level0Id].all_ranks(bitId);
        for (uint64_t symb{0}; symb < TSigma; ++symb) {
            res[symb] =   level0[level0Id/2][symb]
                        + level1[level1Id][symb]
                        + superBlocks[superBlockId][symb]
                        + res[symb]*sign;
            if (symb > 0) {
                prs[symb] = prs[symb-1] + res[symb];
            } else {
                prs[symb] = 0;
            }
        }
        return {res, prs};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
//        ar(bits, level0, level1, superBlocks, totalLength);
        ar(bits, level0, level1, superBlocks, level0_prefix, level1_prefix, superBlocks_prefix, totalLength);

//        std::cout << "bits: " << bits.size() << " " << sizeof(bits[0]) << " -> " << bits.size() * sizeof(bits[0]) << "\n";
//        std::cout << "level0: " << level0.size() << " " << sizeof(level0[0]) << " -> " << level0.size() * sizeof(level0[0]) << "\n";
//        std::cout << "level1: " << level1.size() << " " << sizeof(level1[0]) << " -> " << level1.size() * sizeof(level1[0]) << "\n";
//        std::cout << "superBlocks: " << superBlocks.size() << " " << sizeof(superBlocks[0]) << " -> " << superBlocks.size() * sizeof(superBlocks[0]) << "\n";
    }
};

template <size_t Sigma>
using PartialPairedL1L2L3_64ShortEPRV8 = PartialPairedL1L2L3_NEPRV8<Sigma, 64, uint8_t, uint16_t>;
static_assert(checkRankVector<PartialPairedL1L2L3_64ShortEPRV8>);

template <size_t Sigma>
using PartialPairedL1L2L3_128ShortEPRV8 = PartialPairedL1L2L3_NEPRV8<Sigma, 128, uint8_t, uint16_t>;
static_assert(checkRankVector<PartialPairedL1L2L3_128ShortEPRV8>);

template <size_t Sigma>
using PartialPairedL1L2L3_64EPRV8 = PartialPairedL1L2L3_NEPRV8<Sigma, 64>;
static_assert(checkRankVector<PartialPairedL1L2L3_64EPRV8>);

template <size_t Sigma>
using PartialPairedL1L2L3_128EPRV8 = PartialPairedL1L2L3_NEPRV8<Sigma, 128>;
static_assert(checkRankVector<PartialPairedL1L2L3_128EPRV8>);

template <size_t Sigma>
using PartialPairedL1L2L3_256EPRV8 = PartialPairedL1L2L3_NEPRV8<Sigma, 256>;
static_assert(checkRankVector<PartialPairedL1L2L3_256EPRV8>);

template <size_t Sigma>
using PartialPairedL1L2L3_512EPRV8 = PartialPairedL1L2L3_NEPRV8<Sigma, 512>;
static_assert(checkRankVector<PartialPairedL1L2L3_512EPRV8>);

template <size_t Sigma>
using PartialPairedL1L2L3_1024EPRV8 = PartialPairedL1L2L3_NEPRV8<Sigma, 1024>;
static_assert(checkRankVector<PartialPairedL1L2L3_1024EPRV8>);

template <size_t Sigma>
using PartialPairedL1L2L3_2048EPRV8 = PartialPairedL1L2L3_NEPRV8<Sigma, 2048>;
static_assert(checkRankVector<PartialPairedL1L2L3_2048EPRV8>);


}
