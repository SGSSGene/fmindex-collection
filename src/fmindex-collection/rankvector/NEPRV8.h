// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../bitset_popcount.h"
#include "concepts.h"
#include "utils.h"
#include "EPRV3.h"

#include <bit>
#include <vector>

#if __has_include(<cereal/types/bitset.hpp>)
#include <cereal/types/bitset.hpp>
#endif

#include <iostream>

#include "../ternarylogic.h"

namespace fmindex_collection::rankvector {


namespace neprv8_detail {
    template <size_t N, size_t bitct>
    auto rank(std::array<std::bitset<N>, bitct> const& arr, uint8_t symb) {
        if constexpr (bitct == 3) {
            return mark_exact_fast(symb, arr[2], arr[1], arr[0]);
        } else {
            return mark_exact_large(symb, arr);
        }
    }

    template <size_t All, size_t N, size_t bitct>
    auto rank_all(std::array<std::bitset<N>, bitct> const& arr) {
        if constexpr (bitct == 3) {
            return mark_exact_all(arr[2], arr[1], arr[0]);
        } else {
            auto const& mask = mask_positive_or_negative<N>;
            auto res = std::array<std::bitset<N>, All>{};
            {
                auto p = arr[0] ^ mask[0];
                res[0] = p;
                p = ~p;
                res[1] = p;
            }
/*            p = arr[1] ^ mask[1];
            res[2] = p & res[0];
            res[3] = p & res[1];
            p = ~p;
            res[0] = p & res[0];
            res[1] = p & res[1];

            p = arr[2] ^ mask[1];
            res[4] =  p & res[0];
            res[5] =  p & res[1];
            res[6] =  p & res[2];
            res[7] =  p & res[3];
            p = ~p;
            res[0] = p & res[0];
            res[1] = p & res[1];
            res[2] = p & res[2];
            res[3] = p & res[3];

            p = arr[3] ^ mask[1];
            res[ 8] =  p & res[0];
            res[ 9] =  p & res[1];
            res[10] =  p & res[2];
            res[11] =  p & res[3];
            res[12] =  p & res[4];
            res[13] =  p & res[5];
            res[14] =  p & res[6];
            res[15] =  p & res[7];

            p = ~p;
            res[0] =  p & res[0];
            res[1] =  p & res[1];
            res[2] =  p & res[2];
            res[3] =  p & res[3];
            res[4] =  p & res[4];
            res[5] =  p & res[5];
            res[6] =  p & res[6];
            res[7] =  p & res[7];*/


            auto f = [&](size_t j) {
                auto p = arr[j] ^ mask[1];
                auto range = (1ull<<j);
                for (size_t i{0}; i < range; ++i) { //4-5
                    res[range+i] =  p & res[i];
                }
                p = ~p;
                for (size_t i{0}; i < range; ++i) { // 6-7
                    res[i] = p & res[i];
                }
            };
            for (size_t j{1}; j < bitct; ++j) {
                f(j);
            }
/*            auto res = std::array<std::bitset<N>, All>{};
            for (size_t i{0}; i < All; ++i) {
                res[i] = rank(arr, i);
            }*/
            #if !NDEBUG
            for (size_t i{0}; i < All; ++i) {
                auto v1 = res[i];
                auto v2 = rank(arr, i);
                assert(v1 == v2);
            }
            #endif
            return res;
        }
    }


    template <size_t N, size_t bitct>
    auto prefix_rank(std::array<std::bitset<N>, bitct> const& arr, uint8_t symb) {
        if constexpr (bitct == 3) {
            return mark_exact_or_less_fast(symb, arr[2], arr[1], arr[0]);
        } else {
            return mark_exact_or_less_large(symb, arr);
        }
    }
}
/*
struct FillWidget {
    size_t total{};

    uint64_t currentValue{};

    void addLevel(size_t width) {
    }

    void push(bool v) {
        total += 1;
        currentValue += 
        total += 
    }
};*/

template <size_t TSigma, size_t popcount_width, typename blockL0_t = uint16_t, typename blockL1_t = uint32_t>
struct NEPRV8 {

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
            assert(idx < popcount_width);
            auto v = neprv8_detail::rank(bits, symb);
            return lshift_and_count(v, popcount_width-idx);
        }

        uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
            assert(idx < popcount_width);
            auto v = neprv8_detail::prefix_rank(bits, symb);
            return lshift_and_count(v, popcount_width-idx);
        }

        auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
            assert(idx < popcount_width);

            auto vs = neprv8_detail::rank_all<TSigma>(bits);
            auto v = std::array<uint64_t, TSigma>{};
            static_assert(v.size() <= vs.size());
            for (size_t i{0}; i < v.size(); ++i) {
                v[i] = lshift_and_count(vs[i], popcount_width-idx);
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

    static_assert(popcount_width < (1ull<<level0_size));
    static_assert(level0_size < (1ull<<level1_size));

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

    NEPRV8() = default;
    NEPRV8(std::span<uint8_t const> _symbols) {
        totalLength = _symbols.size();
        auto const length = _symbols.size();
        level1.reserve(length/(1ull<<level1_size)+2);
        bits.reserve(length/popcount_width+2);

        std::array<blockL0_t, TSigma> blockL0_acc{0};
        std::array<blockL1_t, TSigma> blockL1_acc{0};
        std::array<uint64_t, TSigma> sblock_acc{0};

        for (uint64_t size{0}; size < length+1; ++size) {
            if (size % (1ull<<level1_size) == 0) { // new l3 block
                superBlocks.emplace_back(sblock_acc);
                level1.emplace_back();
                bits.emplace_back();
                level0.emplace_back();

                blockL0_acc = {};
                blockL1_acc = {};
            } else if (size % (1ull<<level0_size) == 0) { // new l1 block
                level1.emplace_back(blockL1_acc);
                bits.emplace_back();
                level0.emplace_back();

                blockL0_acc = {};
            } else if (size % popcount_width == 0) { // new l0 block
                bits.emplace_back();
                level0.push_back(blockL0_acc);
            }
            // Abort, we only wanted to add a new block to the end, if required
            if (size == length) continue;

            auto level0Id     = size >>  popcount_width_bits;
            auto bitId        = size &  (popcount_width-1);

            uint64_t symb = _symbols[size];

            for (uint64_t i{}; i < bitct; ++i) {
                auto b = ((symb>>i)&1);
                bits[level0Id].bits[i].set(bitId, b);
            }
            blockL0_acc[symb] += 1;
            blockL1_acc[symb] += 1;
            sblock_acc[symb] += 1;
        }

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
        //prefetch(idx);

        auto level0Id     = idx >> popcount_width_bits;
        auto level1Id     = idx >> level0_size;
        auto superBlockId = idx >> level1_size;
        auto bitId        = idx &  (popcount_width - 1);
        return  size_t{}
                + bits[level0Id].rank(bitId, symb)
                + level0[level0Id][symb]
                + level1[level1Id][symb]
                + superBlocks[superBlockId][symb];
    }

    uint64_t prefix_rank(uint64_t idx, uint8_t symb) const {
        //prefetch(idx);

        auto level0Id     = idx >> popcount_width_bits;
        auto level1Id     = idx >> level0_size;
        auto superBlockId = idx >> level1_size;
        auto bitId        = idx & (popcount_width-1);

        return size_t{}
               + bits[level0Id].prefix_rank(bitId, symb)
               + level0_prefix[level0Id][symb]
               + level1_prefix[level1Id][symb]
               + superBlocks_prefix[superBlockId][symb];

    }


    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
        //prefetch(idx);

        auto level0Id     = idx >> popcount_width_bits;
        auto level1Id     = idx >> level0_size;
        auto superBlockId = idx >> level1_size;
        auto bitId        = idx & (popcount_width-1);

        auto res = bits[level0Id].all_ranks(bitId);
        for (uint64_t symb{0}; symb < TSigma; ++symb) {
            res[symb] =   level0[level0Id][symb]
                        + level1[level1Id][symb]
                        + superBlocks[superBlockId][symb]
                        + res[symb];
        }
        return res;
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, TSigma>, std::array<uint64_t, TSigma>> {
        //prefetch(idx);

        auto rs = all_ranks(idx);
        auto prs = rs;
        for (size_t i{1}; i < prs.size(); ++i) {
            prs[i] = prs[i] + prs[i-1];
        }
        return {rs, prs};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bits, level0, level1, superBlocks, level0_prefix, level1_prefix, superBlocks_prefix, totalLength);
//        ar(bits, level0, level1, superBlocks, totalLength);

//        std::cout << "bits: " << bits.size() << " " << sizeof(bits[0]) << " -> " << bits.size() * sizeof(bits[0]) << "\n";
//        std::cout << "level0: " << level0.size() << " " << sizeof(level0[0]) << " -> " << level0.size() * sizeof(level0[0]) << "\n";
//        std::cout << "level1: " << level1.size() << " " << sizeof(level1[0]) << " -> " << level1.size() * sizeof(level1[0]) << "\n";
//        std::cout << "superBlocks: " << superBlocks.size() << " " << sizeof(superBlocks[0]) << " -> " << superBlocks.size() * sizeof(superBlocks[0]) << "\n";
    }
};

template <size_t Sigma>
using NEPRV8_64 = NEPRV8<Sigma, 64>;
static_assert(checkRankVector<NEPRV8_64>);

template <size_t Sigma>
using NEPRV8_128 = NEPRV8<Sigma, 128>;
static_assert(checkRankVector<NEPRV8_128>);

template <size_t Sigma>
using NEPRV8_256 = NEPRV8<Sigma, 256>;
static_assert(checkRankVector<NEPRV8_256>);

template <size_t Sigma>
using NEPRV8_512 = NEPRV8<Sigma, 512>;
static_assert(checkRankVector<NEPRV8_512>);

template <size_t Sigma>
using NEPRV8_1024 = NEPRV8<Sigma, 1024>;
static_assert(checkRankVector<NEPRV8_1024>);

template <size_t Sigma>
using NEPRV8_2048 = NEPRV8<Sigma, 2048>;
static_assert(checkRankVector<NEPRV8_2048>);


}
