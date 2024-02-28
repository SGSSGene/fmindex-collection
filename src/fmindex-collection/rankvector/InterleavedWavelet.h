// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../builtins.h"
#include "concepts.h"
#include "utils.h"

#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

namespace fmindex_collection::rankvector {

/* Implements the concept `RankVector`
 *
 * \param TSigma size of the alphabet
 */
template <size_t TSigma>
struct InterleavedWavelet {
    static constexpr size_t Sigma = TSigma;

    /** traverse the bits of symb
     */
    template <typename CB>
    static void traverse_symb_bits (uint64_t symb, CB cb) {
        uint64_t id{0};
        uint64_t mask = 1u << (bitct-1);
        for (uint64_t b{0}; b < bitct; ++b) {
            auto bit = (symb & mask) != 0;
            cb(id, bit);
            id = id * 2 + 1 + bit;
            mask = mask >> 1ull;
        }
    }

    // number of full length bitvectors needed `2^bitct ≥ TSigma`
    static constexpr auto bitct = required_bits(TSigma);
    // next full power of 2
    static constexpr auto bvct  = pow(2, bitct);

    struct Block {
        std::array<uint32_t, TSigma> blocks{};
        std::array<uint64_t, bitct>  bits{};

        auto prefetch() const {
            __builtin_prefetch(reinterpret_cast<void const*>(&blocks), 0, 0);
//            __builtin_prefetch((const void*)(&bits), 0, 0);
        }

        uint64_t rank(uint64_t idx, uint64_t symb) const {
            uint64_t s = 0;
            uint64_t l1 = 64;
            uint64_t l2 = idx;
            auto const& bbits = bits;
            uint64_t depth = 0;
            traverse_symb_bits(symb, [&]([[maybe_unused]] uint64_t id, uint64_t bit) {
                auto bits = std::bitset<64>(bbits[depth]) >> s;
                auto c1 = (bits << (64-l1)).count();
                auto d1 = (bits << (64-l2)).count();

                if (bit == 0) {
                    //s = s;
                    l2 = l2 - d1;
                    l1 = l1 - c1;
                } else {
                    s = s + l1 - c1;
                    l2 = d1;
                    l1 = c1;
                }
                depth += 1;
            });
            return l2 + blocks[symb];
        }

        uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
            uint64_t s = 0;
            uint64_t l1 = 64;
            uint64_t l2 = idx;
            uint64_t a{};
            uint64_t depth = 0;
            auto const& bbits = bits;
            traverse_symb_bits(symb+1, [&]([[maybe_unused]] uint64_t id, uint64_t bit) {
                auto bits = std::bitset<64>(bbits[depth]) >> s;
                auto c1 = (bits << (64-l1)).count();
                auto d1 = (bits << (64-l2)).count();

                auto c0 = l1 - c1;
                auto d0 = l2 - d1;
                if (bit == 0) {
                    //s = s;
                    l2 = d0;
                    l1 = c0;
                } else {
                    a += d0;
                    s = s + c0;
                    l2 = d1;
                    l1 = c1;
                }
                depth += 1;
            });
            for (uint64_t i{0}; i <= symb; ++i) {
                a += blocks[i];
            }
            return a;
        }

        template <typename CB>
        void bvAccess(uint64_t idx, CB cb) const {
            bvAccess(0, 64, idx, cb);
        }

        template <uint64_t Depth=0, uint64_t I=0, typename CB>
        void bvAccess(uint64_t s, uint64_t l1, uint64_t idx, CB cb) const {
            auto bv = std::bitset<64>(bits[Depth]) >> s;
            auto d1 = (bv << (64-idx)).count();
            auto d0 = idx - d1;
            if constexpr (Depth < bitct-1) {
                auto c1 = (bv << (64-l1)).count();
                auto c0 = l1 - c1;
                bvAccess<Depth+1, (I << 1)>  (s,    c0, d0, cb);
                bvAccess<Depth+1, (I << 1)+1>(s+c0, c1, d1, cb);
            } else {
                cb(d0, d1, std::integer_sequence<uint64_t, (I << 1)>{});
            }
        }
        template <uint64_t Depth, uint64_t I, typename CB>
        void bvAccess0(CB cb) const {
            if constexpr (Depth < bitct-1) {
                bvAccess0<Depth+1, (I << 1)>(cb);
                bvAccess0<Depth+1, (I << 1)+1>(cb);
            } else {
                cb(0, 0, std::integer_sequence<uint64_t, (I << 1)>{});
            }
        }

        auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
            std::array<uint64_t, TSigma> rs{0};
            bvAccess(idx, [&]<uint64_t I>(uint64_t d0, uint64_t d1, std::integer_sequence<uint64_t, I>) {
                if constexpr (I < TSigma) {
                    rs[I]    = d0 + blocks[I];
                }
                if constexpr (I+1 < TSigma) {
                    rs[I+1]  = d1 + blocks[I+1];
                }
            });
            return rs;
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar(blocks, bits);
        }
    };


    std::vector<Block> blocks;
    std::vector<std::array<uint64_t, TSigma>> superBlocks;
    size_t totalLength{};

    InterleavedWavelet() = default;
    InterleavedWavelet(std::span<uint8_t const> _symbols) {
        totalLength = _symbols.size();
        auto const length = _symbols.size();
        blocks.reserve(length/64+2);

        std::array<uint64_t, bvct> sblock_acc{0};
        std::array<uint32_t, bvct> block_acc{0};
        std::array<uint64_t, bvct> bv_size{};
        for (auto& c : bv_size) {
            c = 0;
        }

        std::array<std::vector<uint8_t>, bvct> waveletcount{};
        auto insertCount = [&]() {
            if (!blocks.empty()) {
                auto& bits = blocks.back().bits;
                uint64_t idx{0};
                for (uint64_t i{0}; i < waveletcount.size(); ++i) {
                    auto t = waveletcount.at(i);
                    for (uint64_t bit : t) {
                        bits[idx / 64] = bits[idx / 64] | (bit << (idx % 64));
                        ++idx;
                    }
                }
            }
            if (blocks.size() % (1ull<<26) == 0) { // new super block + new block
                superBlocks.emplace_back();
                for (uint64_t j{0}; j < TSigma; ++j) {
                    superBlocks.back()[j] = sblock_acc[j];
                }
                block_acc = {};
            }

            blocks.emplace_back();
            for (uint64_t i{0}; i < TSigma; ++i) {
                blocks.back().blocks[i] = block_acc[i];
            }
            for (auto& c : waveletcount) {
                c.clear();
            }
        };

        for (uint64_t size{0}; size < length; ++size) {
            if (size % 64 == 0) { // new block
               insertCount();
            }

            auto symb = _symbols[size];
            traverse_symb_bits(symb, [&](uint64_t id, uint64_t bit) {
                waveletcount[id].push_back(bit);
            });
            block_acc[symb] += 1;
            sblock_acc[symb] += 1;
        }
        while (waveletcount[0].size() < 64) {
            traverse_symb_bits(bvct-1, [&](uint64_t id, uint64_t bit) {
                waveletcount[id].push_back(bit);
            });
            block_acc[bvct-1] += 1;
            sblock_acc[bvct-1] += 1;
        }
        insertCount();
    }

    size_t size() const {
        return totalLength;
    }

    auto prefetch(uint64_t idx) const {
        auto blockId      = idx >>  6;
        blocks[blockId].prefetch();
    }

    uint8_t symbol(uint64_t idx) const {
        //!TODO needs a better implementation than this....
        for (size_t i{0}; i < Sigma; ++i) {
            if (rank(idx, i) != rank(idx+1, i)) {
                return i;
            }
        }
        throw std::runtime_error{"this rank vector seems to be in an invalid state"};
    }

    uint64_t rank(uint64_t idx, uint8_t symb) const {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> 32;
        auto bitId        = idx &  63;
        return blocks[blockId].rank(bitId, symb) + superBlocks[superBlockId][symb];
    }

    uint64_t prefix_rank(uint64_t idx, uint8_t symb) const {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> 32;
        auto bitId        = idx &  63;
        uint64_t a={};
        for (uint64_t i{0}; i<= symb; ++i) {
            a += superBlocks[superBlockId][i];
        }
        return blocks[blockId].prefix_rank(bitId, symb) + a;
    }

    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> 32;
        auto bitId        = idx &  63;

        auto rs = blocks[blockId].all_ranks(bitId);

        rs[0] += superBlocks[superBlockId][0];
        for (uint64_t i{1}; i < TSigma; ++i) {
            rs[i] += superBlocks[superBlockId][i];
        }

        return rs;
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, TSigma>, std::array<uint64_t, TSigma>> {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> 32;
        auto bitId        = idx &  63;

        auto rs = blocks[blockId].all_ranks(bitId);

        std::array<uint64_t, TSigma> prs;

        rs[0] += superBlocks[superBlockId][0];
        prs[0] = rs[0];
        for (uint64_t i{1}; i < TSigma; ++i) {

            rs[i] += superBlocks[superBlockId][i];
            prs[i] = prs[i-1] + rs[i];
        }

        return {rs, prs};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(blocks, superBlocks, totalLength);
    }
};

static_assert(checkRankVector<InterleavedWavelet>);

}
