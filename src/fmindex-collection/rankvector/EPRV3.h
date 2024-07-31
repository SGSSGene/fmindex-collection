// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "concepts.h"
#include "utils.h"

#include <bitset>
#include <cassert>
#include <vector>

namespace fmindex_collection::rankvector {

template <size_t TSigma, typename block_t = uint16_t>
struct EPRV3 {

    static constexpr size_t Sigma = TSigma;

    // number of full length bitvectors needed `2^bitct ≥ TSigma`
    static constexpr auto bitct = required_bits(TSigma-1);
    // next full power of 2
    static constexpr auto bvct  = pow(2, bitct);

    struct Block {
        std::array<block_t, TSigma> blocks{};

        auto operator[](uint8_t symb) const {
            return blocks[symb];
        }
        auto& operator[](uint8_t symb) {
            return blocks[symb];
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar(blocks);
        }
    };

    struct InBits {
        std::array<uint64_t, bitct> bits{};

        uint8_t symbol(uint64_t idx) const {
            uint8_t symb{};
            for (uint64_t i{bitct}; i > 0; --i) {
                auto b = (bits[i-1] >> idx) & 1;
                symb = (symb<<1) | b;
            }
            return symb;
        }

        uint64_t rank(uint64_t idx, uint8_t symb) const {
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

            auto rs = std::array<uint64_t, TSigma>{0};

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

    static constexpr uint64_t block_size = sizeof(block_t) * 8;

    std::vector<InBits> bits;
    std::vector<Block> blocks_;
    std::vector<std::array<uint64_t, TSigma>> superBlocks;
    size_t totalLength;

    EPRV3() = default;
    EPRV3(std::span<uint8_t const> _symbols) {
        totalLength = _symbols.size();
        blocks_.reserve(_symbols.size()/block_size+1);
        bits.reserve(_symbols.size()/block_size+1);

        auto sblock_acc = std::array<uint64_t, TSigma>{}; // accumulator for super blocks
        auto block_acc  = std::array<block_t, TSigma>{};  // accumulator for blocks

        for (uint64_t size{0}; size < _symbols.size();) {
            superBlocks.emplace_back(sblock_acc);
            block_acc = {};

             for (uint64_t blockId{0}; blockId < (uint64_t{1}<<block_size)/64 and size < _symbols.size(); ++blockId) {
                blocks_.emplace_back();
                blocks_.back().blocks = block_acc;
                bits.emplace_back();

                for (uint64_t bitId{0}; bitId < 64 and size < _symbols.size(); ++bitId, ++size) {

                    uint64_t symb = _symbols[size];

                    for (uint64_t i{}; i < bitct; ++i) {
                        auto b = ((symb>>i)&1);
                        bits.back().bits[i] |= (b << bitId);
                    }

                    block_acc[symb] += 1;
                    sblock_acc[symb] += 1;
                }
            }
        }
        // For safety we add a new super block and block
        superBlocks.emplace_back(sblock_acc);
        blocks_.emplace_back();
        blocks_.back().blocks = block_acc;
        bits.emplace_back();
    }

    void prefetch(uint64_t idx) const {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;

        __builtin_prefetch(reinterpret_cast<void const*>(&blocks_[blockId]), 0, 0);
        __builtin_prefetch(reinterpret_cast<void const*>(&bits[blockId]), 0, 0);
        __builtin_prefetch(reinterpret_cast<void const*>(&superBlocks[superBlockId]), 0, 0);
//        blocks[blockId].prefetch();
    }

    size_t size() const {
        return totalLength;
    }

    uint8_t symbol(uint64_t idx) const {
        auto blockId      = idx >>  6;
        auto bitId        = idx &  63;
        return bits[blockId].symbol(bitId);
    }

    uint64_t rank(uint64_t idx, uint8_t symb) const {
        prefetch(idx);

        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;
        auto bitId        = idx &  63;
        return blocks_[blockId][symb] + bits[blockId].rank(bitId, symb) + superBlocks[superBlockId][symb];
    }

    uint64_t prefix_rank(uint64_t idx, uint8_t symb) const {
        prefetch(idx);

        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;
        auto bitId        = idx &  63;
        uint64_t a={};
        for (uint64_t i{0}; i<= symb; ++i) {
            a += superBlocks[superBlockId][i] + blocks_[blockId][i];
        }
        return bits[blockId].prefix_rank(bitId, symb) + a;
    }


    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
        prefetch(idx);

        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;
        auto bitId        = idx &  63;
        auto res = std::array<uint64_t, TSigma>{};
        for (uint64_t symb{0}; symb < TSigma; ++symb) {
            res[symb] = blocks_[blockId][symb] + bits[blockId].rank(bitId, symb) + superBlocks[superBlockId][symb];
        }
        return res;
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, TSigma>, std::array<uint64_t, TSigma>> {
        prefetch(idx);

        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;
        auto bitId        = idx &  63;

        std::array<uint64_t, TSigma> prs;
        auto rs = bits[blockId].all_ranks(bitId);

        rs[0] += superBlocks[superBlockId][0] + blocks_[blockId][0];
        prs[0] = rs[0];
        for (uint64_t symb{1}; symb < TSigma; ++symb) {
            prs[symb] = prs[symb-1] + superBlocks[superBlockId][symb] + rs[symb] + blocks_[blockId][symb];
            rs[symb] += superBlocks[superBlockId][symb] + blocks_[blockId][symb];
        }
        return {rs, prs};
    }

    uint64_t rank_symbol(uint64_t idx) const {
        prefetch(idx);

        auto blockId = idx >> 6;
        auto bitId = idx & 63;
        auto superBlockId = idx >> block_size;
        auto [rank, symb] = bits[blockId].rank_symbol(bitId);
        return rank + superBlocks[superBlockId][symb] + blocks_[blockId][symb];
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(blocks_, bits, superBlocks, totalLength);
    }

};

template <size_t TSigma> using EPRV3_8  = EPRV3<TSigma,  uint8_t>;
static_assert(checkRankVector<EPRV3_8>);

template <size_t TSigma> using EPRV3_16 = EPRV3<TSigma, uint16_t>;
static_assert(checkRankVector<EPRV3_16>);

#if SIZE_MAX == UINT64_MAX
template <size_t TSigma> using EPRV3_32 = EPRV3<TSigma, uint32_t>;
static_assert(checkRankVector<EPRV3_32>);
#endif

}
