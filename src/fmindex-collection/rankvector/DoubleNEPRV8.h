// SPDX-FileCopyrightText: 2024  Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "concepts.h"
#include "utils.h"
#include "EPRV3.h"

#include <vector>

#if __has_include(<cereal/types/bitset.hpp>)
#include <cereal/types/bitset.hpp>
#endif

#include <iostream>

namespace fmindex_collection::rankvector {

template <size_t TSigma, size_t popcount_width, typename blockL0_t = uint16_t, typename blockL1_t = uint32_t>
struct DoubleNEPRV8 {

    static constexpr size_t Sigma = TSigma;

    // number of full length bitvectors needed `2^bitct ≥ TSigma`
    static constexpr auto bitct = required_bits(TSigma-1);
    // next full power of 2
    static constexpr auto bvct  = pow(2, bitct);

    static constexpr auto popcount_width_bits = required_bits(popcount_width-1);

    struct InBits {
        std::array<std::bitset<popcount_width>, bitct> bits;

        uint8_t symbol(uint64_t idx) const {
            uint8_t symb{};
            for (uint64_t i{bitct}; i > 0; --i) {
                auto b = (bits[i-1] >> idx).test(0);
                symb = (symb<<1) | b;
            }
            return symb;
        }

        template <bool reverse=false>
        uint64_t rank(uint64_t idx, uint8_t symb) const {
            assert(idx < popcount_width);
            auto f = [&]<uint64_t I>(std::integer_sequence<uint64_t, I>) {
                if (symb & (1<<I)) {
                    return bits[I];
                } else {
                    return ~bits[I];
                }
            };
            auto mask = [&]<uint64_t ...Is>(std::integer_sequence<uint64_t, Is...>) {
                return (f(std::integer_sequence<uint64_t, Is>{})&...);
            }(std::make_integer_sequence<uint64_t, bitct>{});

            if constexpr (reverse) {
                auto bitset = mask >> idx;
                return bitset.count();
            }
            auto bitset = mask << (popcount_width-idx);
            return bitset.count();
        }

        template <bool reverse=false>
        uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
            assert(idx < popcount_width);
            auto f = [&]<uint64_t I>(std::integer_sequence<uint64_t, I>, uint64_t _symb) {
                if (_symb & (1<<I)) {
                    return bits[I];
                } else {
                    return ~bits[I];
                }
            };
            auto mask = std::bitset<popcount_width>{};

            for (uint64_t i{0}; i <= symb; ++i) {
                mask |= [&]<uint64_t ...Is>(std::integer_sequence<uint64_t, Is...>) {
                    return (f(std::integer_sequence<uint64_t, Is>{}, i)&...);
                }(std::make_integer_sequence<uint64_t, bitct>{});
            }

            if constexpr (reverse) {
                auto bitset = mask >> idx;
                return bitset.count();
            }

            auto bitset = mask << (popcount_width-idx);
            return bitset.count();
        }

        template <bool reverse=false>
        auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
            assert(idx < popcount_width);

            auto rs = std::array<uint64_t, TSigma>{0};

            auto f = [&]<uint64_t I>(uint64_t _symb, std::integer_sequence<uint64_t, I>) {
                if (_symb & (1<<I)) {
                    return bits[I];
                } else {
                    return ~bits[I];
                }
            };

            for (uint64_t i{0}; i < TSigma; ++i) {
                auto mask = [&]<uint64_t ...Is>(std::integer_sequence<uint64_t, Is...>) {
                    return (f(i, std::integer_sequence<uint64_t, Is>{})&...);
                }(std::make_integer_sequence<uint64_t, bitct>{});
                if constexpr (reverse) {
                    rs[i] = (mask >> idx).count();
                } else {
                    rs[i] = (mask << (popcount_width - idx)).count();
                }
            }
            return rs;
        }

        template <bool reverse=false>
        auto rank_symbol(uint64_t idx) const -> std::tuple<uint64_t, uint64_t> {
            assert(idx < popcount_width);

            uint64_t symb{};
            auto mask = std::bitset<popcount_width>{};
            auto f = [&]<uint64_t I>(std::integer_sequence<uint64_t, I>) {
                auto b = (bits[I] >> idx).test(0);
                mask |= bits[I] ^ -b;
                symb |= b << I;
            };
            [&]<uint64_t ...Is>(std::integer_sequence<uint64_t, Is...>) {
                (f(std::integer_sequence<uint64_t, Is>{}) ,...);
            }(std::make_integer_sequence<uint64_t, bitct>{});

            if constexpr (reverse) {
                auto bitset = (~mask) >> idx;
                return {bitset.count(), symb};
            }

            auto bitset = (~mask) << (popcount_width-idx);
            return {bitset.count(), symb};
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar(bits);
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

    size_t totalLength;

    DoubleNEPRV8() = default;
    DoubleNEPRV8(std::span<uint8_t const> _symbols) {
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

/*            {
                size_t a{};
                for (auto v : level0.back()) {
                    a += v;
                }
                assert(a == 64+(((level0.size()-1)%2)*128));
            }*/
        }

/*        for (size_t i{0}; i < level0.size(); ++i) {
            size_t a{};
            for (auto v : level0[i]) {
                a += v;
            }
            assert(a == 64+((i%2)*128));
        }*/

        // check last level0, this should be more efficient
/*        {
            uint64_t a{};
            for (size_t symb{}; symb < Sigma; ++symb) {
                a += level0.back()[symb];
            }
            assert(a <= popcount_width);
            level0.back()[0] += (popcount_width - a);
        }*/
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
        auto bitId        = idx &  (popcount_width - 1);
        if ((idx / popcount_width) % 2 == 1) {
            return  size_t{}
                    + bits[level0Id].rank(bitId, symb)
                    + level0[level0Id/2][symb]
                    + level1[level1Id][symb]
                    + superBlocks[superBlockId][symb];
        } else {
            return  size_t{}
                    + level0[level0Id/2][symb]
                    + level1[level1Id][symb]
                    + superBlocks[superBlockId][symb]
                    - bits[level0Id].template rank</*reverse = */true>(bitId, symb);
        }
    }

    uint64_t prefix_rank(uint64_t idx, uint8_t symb) const {
        prefetch(idx);

        auto level0Id     = idx >>  popcount_width_bits;
        auto level1Id     = idx >> level0_size;
        auto superBlockId = idx >> level1_size;
        auto bitId        = idx &  (popcount_width - 1);
        uint64_t a={};

        if ((idx / popcount_width) % 2 == 1) {
            for (uint64_t i{0}; i <= symb; ++i) {
                a += size_t{}
                     + level0[level0Id/2][i]
                     + level1[level1Id][i]
                     + superBlocks[superBlockId][i];

            }
            return a + bits[level0Id].prefix_rank(bitId, symb);
        } else {
            for (uint64_t i{0}; i <= symb; ++i) {
                a += size_t{}
                     + level0[level0Id/2][i]
                     + level1[level1Id][i]
                     + superBlocks[superBlockId][i];

            }
            auto c = bits[level0Id].template prefix_rank</*reverse = */true>(bitId, symb);
            return a - c;

        }
    }


    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
        prefetch(idx);

        auto level0Id     = idx >>  popcount_width_bits;
        auto level1Id     = idx >> level0_size;
        auto superBlockId = idx >> level1_size;
        auto bitId        = idx &  (popcount_width - 1);
        auto res = std::array<uint64_t, TSigma>{};
        if ((idx / popcount_width) % 2 == 1) {
            for (uint64_t symb{0}; symb < TSigma; ++symb) {
                res[symb] =   bits[level0Id].rank(bitId, symb)
                            + level0[level0Id/2][symb]
                            + level1[level1Id][symb]
                            + superBlocks[superBlockId][symb];
            }
        } else {
            for (uint64_t symb{0}; symb < TSigma; ++symb) {
                res[symb] =   level0[level0Id/2][symb]
                            + level1[level1Id][symb]
                            + superBlocks[superBlockId][symb]
                            - bits[level0Id].template rank</*reverse =*/ true>(bitId, symb);
            }

        }
        return res;
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, TSigma>, std::array<uint64_t, TSigma>> {
        prefetch(idx);

        auto level0Id     = idx >>  popcount_width_bits;
        auto level1Id     = idx >> level0_size;
        auto superBlockId = idx >> level1_size;
        auto bitId        = idx &  (popcount_width - 1);

        auto prs = std::array<uint64_t, TSigma>{};

        if ((idx / popcount_width) % 2 >= 1) {
            auto rs = bits[level0Id].all_ranks(bitId);

            rs[0] +=   level0[level0Id/2][0]
                     + level1[level1Id][0]
                     + superBlocks[superBlockId][0];

            prs[0] = rs[0];
            for (uint64_t symb{1}; symb < TSigma; ++symb) {
                auto a =   level0[level0Id/2][symb]
                         + level1[level1Id][symb]
                         + superBlocks[superBlockId][symb]
                         + rs[symb];

                rs[symb]  = a;
                prs[symb] = prs[symb-1] + a;
            }
            return {rs, prs};
        } else {
            auto rs = bits[level0Id].template all_ranks</*reverse =*/true>(bitId);

            rs[0]  =   level0[level0Id/2][0]
                     + level1[level1Id][0]
                     + superBlocks[superBlockId][0]
                     - rs[0];

            prs[0] = rs[0];
            for (uint64_t symb{1}; symb < TSigma; ++symb) {
                auto a =   level0[level0Id/2][symb]
                         + level1[level1Id][symb]
                         + superBlocks[superBlockId][symb]
                         - rs[symb];

                rs[symb] = a;
                prs[symb] = prs[symb-1] + a;
            }
            return {rs, prs};
        }
    }

    uint64_t rank_symbol(uint64_t idx) const {
        prefetch(idx);

        auto level0Id     = idx >>  popcount_width_bits;
        auto level1Id     = idx >> level0_size;
        auto superBlockId = idx >> level1_size;
        auto bitId        = idx & (popcount_width-1);

        if ((idx / popcount_width) % 2 == 1) {
            auto [rank, symb] = bits[level0Id].rank_symbol(bitId);
            return    rank
                    + level0[level0Id/2][symb]
                    + level1[level1Id][symb]
                    + superBlocks[superBlockId][symb];
        } else {
            auto [rank, symb] = bits[level0Id].template rank_symbol</*reverse = */true>(bitId);
            return    level0[level0Id/2][symb]
                    + level1[level1Id][symb]
                    + superBlocks[superBlockId][symb];
                    - rank;

        }
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bits, level0, level1, superBlocks, totalLength);
//        std::cout << "bits: " << bits.size() << " " << sizeof(bits[0]) << " -> " << bits.size() * sizeof(bits[0]) << "\n";
//        std::cout << "level0: " << level0.size() << " " << sizeof(level0[0]) << " -> " << level0.size() * sizeof(level0[0]) << "\n";
//        std::cout << "level1: " << level1.size() << " " << sizeof(level1[0]) << " -> " << level1.size() * sizeof(level1[0]) << "\n";
//        std::cout << "superBlocks: " << superBlocks.size() << " " << sizeof(superBlocks[0]) << " -> " << superBlocks.size() * sizeof(superBlocks[0]) << "\n";
    }
};

template <size_t Sigma>
using Double64ShortEPRV8 = DoubleNEPRV8<Sigma, 64, uint8_t, uint16_t>;
static_assert(checkRankVector<Double64ShortEPRV8>);

template <size_t Sigma>
using Double128ShortEPRV8 = DoubleNEPRV8<Sigma, 128, uint8_t, uint16_t>;
static_assert(checkRankVector<Double128ShortEPRV8>);

template <size_t Sigma>
using Double64EPRV8 = DoubleNEPRV8<Sigma, 64>;
static_assert(checkRankVector<Double64EPRV8>);

template <size_t Sigma>
using Double128EPRV8 = DoubleNEPRV8<Sigma, 128>;
static_assert(checkRankVector<Double128EPRV8>);

template <size_t Sigma>
using Double256EPRV8 = DoubleNEPRV8<Sigma, 256>;
static_assert(checkRankVector<Double256EPRV8>);

template <size_t Sigma>
using Double512EPRV8 = DoubleNEPRV8<Sigma, 512>;
static_assert(checkRankVector<Double512EPRV8>);

}
