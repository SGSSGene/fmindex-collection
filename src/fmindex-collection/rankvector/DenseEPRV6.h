// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../DenseVector.h"
#include "concepts.h"
#include "utils.h"

#include <bitset>
#include <vector>

namespace fmindex_collection::rankvector {

template <size_t TSigma>
struct DenseEPRV6 {

    // number of full length bitvectors needed `2^bitct ≥ TSigma`
    static constexpr auto bitct = required_bits(TSigma-1);
    // next full power of 2
    static constexpr auto bvct  = pow(2, bitct);

    struct InBits {
        std::array<uint64_t, bitct> bits{};

        uint64_t rank(uint64_t idx, uint64_t symb) const {
            assert(idx < 64);
            auto f = [&]<uint64_t I>(std::integer_sequence<uint64_t, I>) {
                return bits[I] ^ -((~symb>>I)&1);
            };
            auto mask = [&]<uint64_t ...Is>(std::integer_sequence<uint64_t, Is...>) {
                return (f(std::integer_sequence<uint64_t, Is>{})&...);
            }(std::make_integer_sequence<uint64_t, bitct>{});

            auto bitset = std::bitset<64>(mask) << (64-idx);
            return bitset.count();
        }

        uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
            auto f = [&]<uint64_t I>(std::integer_sequence<uint64_t, I>, uint64_t _symb) {
                return bits[I] ^ -((~_symb>>I)&1);
            };
            uint64_t mask{};

            for (uint64_t i{0}; i <= symb; ++i) {
                mask |= [&]<uint64_t ...Is>(std::integer_sequence<uint64_t, Is...>) {
                    return (f(std::integer_sequence<uint64_t, Is>{}, i)&...);
                }(std::make_integer_sequence<uint64_t, bitct>{});
            }

            auto bitset = std::bitset<64>(mask) << (64-idx);
            return bitset.count();
        }

        auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
            assert(idx < 64);

            std::array<uint64_t, TSigma> rs{0};

            auto f = [&]<uint64_t I>(uint64_t symb, std::integer_sequence<uint64_t, I>) {
                return bits[I] ^ -((~symb>>I)&1);
            };

            for (uint64_t i{0}; i < TSigma; ++i) {
                auto mask = [&]<uint64_t ...Is>(std::integer_sequence<uint64_t, Is...>) {
                    return (f(i, std::integer_sequence<uint64_t, Is>{})&...);
                }(std::make_integer_sequence<uint64_t, bitct>{});
                rs[i] = (std::bitset<64>(mask) << (64 - idx)).count();
            }
            return rs;
        }

        uint8_t symbol(uint64_t idx) const {
            uint8_t symb{};
            for (uint64_t i{bitct}; i > 0; --i) {
                auto b = (bits[i-1] >> idx) & 1;
                symb = (symb<<1) | b;
            }
            return symb;
        }

        auto rank_symbol(uint64_t idx) const -> std::tuple<uint64_t, uint64_t> {
            assert(idx < 64);

            uint64_t symb{};
            uint64_t mask{};
            auto f = [&]<uint64_t I>(std::integer_sequence<uint64_t, I>) {
                auto b = (bits[I] >> idx) & 1;
                mask |= bits[I] ^ -b;
                symb |= b << I;
            };
            [&]<uint64_t ...Is>(std::integer_sequence<uint64_t, Is...>) {
                (f(std::integer_sequence<uint64_t, Is>{}) ,...);
            }(std::make_integer_sequence<uint64_t, bitct>{});

            auto bitset = std::bitset<64>{~mask} << (64-idx);

            return {bitset.count(), symb};
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar(bits);
        }
    };

    static constexpr size_t Sigma = TSigma;

    using blockL0_t = uint8_t;
    using blockL1_t = uint16_t;
    static constexpr uint64_t level0_size = sizeof(blockL0_t) * 8;
    static constexpr uint64_t level1_size = sizeof(blockL1_t) * 8;

    using BlockL0 = std::array<blockL0_t, TSigma>;
    using BlockL1 = std::array<blockL1_t, TSigma>;

    std::vector<InBits> bits;
    std::vector<BlockL0> level0;
    std::vector<BlockL1> level1;
    DenseVector superBlocks;

    size_t totalSize;

    DenseEPRV6() = default;
    DenseEPRV6(std::span<uint8_t const> _symbols)
        : superBlocks(std::max(uint64_t{1}, uint64_t(std::ceil(std::log2(_symbols.size())))))
    {
        auto const length = _symbols.size();
        totalSize = length;
        level1.reserve(length/(1ull<<level1_size)+2);
        level0.reserve(length/64+2);
        bits.reserve(length/64+2);

        std::array<blockL0_t, TSigma> blockL0_acc{0};
        std::array<blockL1_t, TSigma> blockL1_acc{0};
        std::array<uint64_t, TSigma> sblock_acc{0};


        for (uint64_t size{0}; size < length+1; ++size) {
            if (size % (1ull<<level1_size) == 0) { // new l3 block
                for (auto v : sblock_acc) {
                    superBlocks.push_back(v);
                }
                level1.emplace_back();
                level0.emplace_back();
                bits.emplace_back();
                blockL0_acc = {};
                blockL1_acc = {};
            } else if (size % (1ull<<level0_size) == 0) { // new l1 block
                level1.emplace_back(blockL1_acc);
                level0.emplace_back();
                bits.emplace_back();

                blockL0_acc = {};
            } else if (size % 64 == 0) { // new l0 block
                level0.emplace_back(blockL0_acc);
                bits.emplace_back();
            }
            // Abort, we only wanted to add a new block to the end, if required
            if (size == length) continue;

            auto level0Id     = size >>  6;
            auto bitId        = size &  63;

            uint64_t symb = _symbols[size];

            for (uint64_t i{}; i < bitct; ++i) {
                auto b = ((symb>>i)&1);
                bits[level0Id].bits[i] |= (b << bitId);
            }
            blockL0_acc[symb] += 1;
            blockL1_acc[symb] += 1;
            sblock_acc[symb] += 1;
        }
    }

    void prefetch(uint64_t idx) const {
        auto level0Id     = idx >>  6;
        auto level1Id     = idx >> level0_size;

        __builtin_prefetch(reinterpret_cast<void const*>(&level1[level1Id]), 0, 0);
        __builtin_prefetch(reinterpret_cast<void const*>(&level0[level0Id]), 0, 0);
        __builtin_prefetch(reinterpret_cast<void const*>(&bits[level1Id]), 0, 0);
    }

    size_t size() const {
        return totalSize;
    }

    uint8_t symbol(uint64_t idx) const {
        auto level0Id     = idx >>  6;
        auto bitId        = idx &  63;
        return bits[level0Id].symbol(bitId);
    }

    uint64_t rank(uint64_t idx, uint8_t symb) const {
        prefetch(idx);

        auto level0Id     = idx >>  6;
        auto level1Id     = idx >> level0_size;
        auto superBlockId = idx >> level1_size;
        auto bitId        = idx &  63;
        return    bits[level0Id].rank(bitId, symb)
                + level0[level0Id][symb]
                + level1[level1Id][symb]
                + superBlocks[superBlockId*TSigma+symb];
    }

    uint64_t prefix_rank(uint64_t idx, uint8_t symb) const {
        prefetch(idx);

        auto level0Id     = idx >>  6;
        auto level1Id     = idx >> level0_size;
        auto superBlockId = idx >> level1_size;
        auto bitId        = idx &  63;
        uint64_t a={};
        for (uint64_t i{0}; i<= symb; ++i) {
            a +=   level0[level0Id][i]
                 + level1[level1Id][i]
                 + superBlocks[superBlockId*TSigma+i];
        }
        return bits[level0Id].prefix_rank(bitId, symb) + a;
    }

    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
        prefetch(idx);

        auto level0Id     = idx >>  6;
        auto level1Id     = idx >> level0_size;
        auto superBlockId = idx >> level1_size;
        auto bitId        = idx &  63;
        auto res = std::array<uint64_t, TSigma>{};
        for (uint64_t symb{0}; symb < TSigma; ++symb) {
            res[symb] =   bits[level0Id].rank(bitId, symb)
                        + level0[level0Id][symb]
                        + level1[level1Id][symb]
                        + superBlocks[superBlockId*TSigma+symb];
        }
        return res;
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, TSigma>, std::array<uint64_t, TSigma>> {
        prefetch(idx);

        auto level0Id     = idx >>  6;
        auto level1Id     = idx >> level0_size;
        auto superBlockId = idx >> level1_size;
        auto bitId        = idx &  63;

        std::array<uint64_t, TSigma> prs;
        auto rs = bits[level0Id].all_ranks(bitId);

        rs[0] +=   level0[level0Id][0]
                 + level1[level1Id][0]
                 + superBlocks[superBlockId*TSigma+0];

        prs[0] = rs[0];
        for (uint64_t symb{1}; symb < TSigma; ++symb) {
            auto a =   level0[level0Id][symb]
                     + level1[level1Id][symb]
                     + superBlocks[superBlockId*TSigma+symb];

            prs[symb] = prs[symb-1] + rs[symb] + a;
            rs[symb] += a;
        }
        return {rs, prs};
    }

    uint64_t rank_symbol(uint64_t idx) const {
        prefetch(idx);

        auto level0Id     = idx >>  6;
        auto level1Id     = idx >> level0_size;
        auto superBlockId = idx >> level1_size;
        auto bitId = idx & 63;

        auto [rank, symb] = bits[level0Id].rank_symbol(bitId);
        return    rank
                + level0[level0Id][symb]
                + level1[level1Id][symb]
                + superBlocks[superBlockId*TSigma+symb];
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bits, level0, level1, superBlocks, totalSize);
    }
};

static_assert(checkRankVector<DenseEPRV6>);

}
