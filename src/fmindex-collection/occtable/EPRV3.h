// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
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
#include <vector>


namespace fmindex_collection {
namespace occtable {
namespace eprV3_impl {

template <uint64_t TSigma, typename block_t>
struct Bitvector {

    // number of full length bitvectors needed `2^bitct ≥ TSigma`
    static constexpr auto bitct = required_bits(TSigma-1);
    // next full power of 2
    static constexpr auto bvct  = pow(2, bitct);

    struct Block {
        std::array<block_t, TSigma> blocks{};

        auto operator[](uint64_t symb) const {
            return blocks[symb];
        }
        auto& operator[](uint64_t symb) {
            return blocks[symb];
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar(blocks);
        }
    };

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

        uint64_t symbol(uint64_t idx) const {
            uint64_t symb{};
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

    static constexpr uint64_t block_size = sizeof(block_t) * 8;

    std::vector<InBits> bits;
    std::vector<Block> blocks_;
    std::vector<std::array<uint64_t, TSigma>> superBlocks;
    std::array<uint64_t, TSigma+1> C;

    Bitvector(std::span<uint8_t const> _bwt) {
        blocks_.reserve(_bwt.size()/block_size+1);
        bits.reserve(_bwt.size()/block_size+1);

        auto sblock_acc = std::array<uint64_t, TSigma>{}; // accumulator for super blocks
        auto block_acc  = std::array<block_t, TSigma>{};  // accumulator for blocks

        for (uint64_t size{0}; size < _bwt.size();) {
            superBlocks.emplace_back(sblock_acc);
            block_acc = {};

             for (uint64_t blockId{0}; blockId < (1ull<<block_size)/64 and size < _bwt.size(); ++blockId) {
                blocks_.emplace_back();
                blocks_.back().blocks = block_acc;
                bits.emplace_back();

                for (uint64_t bitId{0}; bitId < 64 and size < _bwt.size(); ++bitId, ++size) {

                    uint64_t symb = _bwt[size];

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

        C[0] = 0;
        for (uint64_t i{0}; i < TSigma; ++i) {
            C[i+1] = sblock_acc[i] + C[i];
        }
    }

    Bitvector(cereal_tag) {}


    uint64_t memoryUsage() const {
        return blocks_.size() * sizeof(blocks_.back())
            + bits.size() * sizeof(bits)
            + superBlocks.size() * sizeof(superBlocks.back())
            + sizeof(C);
    }

    void prefetch(uint64_t idx) const {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;

        __builtin_prefetch(reinterpret_cast<void const*>(&blocks_[blockId]), 0, 0);
        __builtin_prefetch(reinterpret_cast<void const*>(&bits[blockId]), 0, 0);
        __builtin_prefetch(reinterpret_cast<void const*>(&superBlocks[superBlockId]), 0, 0);
//        blocks[blockId].prefetch();
    }

    uint64_t rank(uint64_t idx, uint64_t symb) const {
        prefetch(idx);

        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;
        auto bitId        = idx &  63;
        return blocks_[blockId][symb] + bits[blockId].rank(bitId, symb) + superBlocks[superBlockId][symb] + C[symb];
    }

    uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
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
            res[symb] = blocks_[blockId][symb] + bits[blockId].rank(bitId, symb) + superBlocks[superBlockId][symb] + C[symb];
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

        rs[0] += superBlocks[superBlockId][0] + C[0] + blocks_[blockId][0];
        prs[0] = rs[0];
        for (uint64_t symb{1}; symb < TSigma; ++symb) {
            prs[symb] = prs[symb-1] + superBlocks[superBlockId][symb] + rs[symb] + blocks_[blockId][symb];
            rs[symb] += C[symb] + superBlocks[superBlockId][symb] + blocks_[blockId][symb];
        }
        return {rs, prs};
    }

    uint64_t symbol(uint64_t idx) const {
        auto blockId      = idx >>  6;
        auto bitId        = idx &  63;
        return bits[blockId].symbol(bitId);
    }

    uint64_t rank_symbol(uint64_t idx) const {
        prefetch(idx);

        auto blockId = idx >> 6;
        auto bitId = idx & 63;
        auto superBlockId = idx >> block_size;
        auto [rank, symb] = bits[blockId].rank_symbol(bitId);
        return rank + superBlocks[superBlockId][symb] + C[symb] + blocks_[blockId][symb];
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(blocks_, bits, superBlocks, C);
    }
};


template <uint64_t TSigma, typename block_t>
struct OccTable {
    using TLengthType = uint64_t;
    static constexpr uint64_t Sigma = TSigma;

    Bitvector<Sigma, block_t> bitvector;

    static uint64_t expectedMemoryUsage(uint64_t length) {
        using Block = typename Bitvector<TSigma, block_t>::Block;
        auto blockSize = std::max(alignof(Block), sizeof(Block));

        uint64_t C           = sizeof(uint64_t) * (Sigma+1);
        uint64_t blocks      = blockSize        * (length+1) / 64;
        uint64_t superblocks = sizeof(uint64_t) * (length+1) / (1ull << (sizeof(block_t) * 8));
        return C + blocks + superblocks;
    }

    OccTable(std::span<uint8_t const> _bwt)
        : bitvector{_bwt}
    {}

    OccTable(cereal_tag)
        : bitvector{cereal_tag{}}
    {}

    uint64_t memoryUsage() const {
        return bitvector.memoryUsage() + sizeof(OccTable);
    }

    uint64_t size() const {
        return bitvector.C.back();
    }

    auto prefetch(uint64_t idx) const {
        bitvector.prefetch(idx);
    }

    uint64_t rank(uint64_t idx, uint64_t symb) const {
        return bitvector.rank(idx, symb);
    }

    uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
        return bitvector.prefix_rank(idx, symb);
    }

    uint64_t symbol(uint64_t idx) const {
        return bitvector.symbol(idx);
    }

    uint64_t rank_symbol(uint64_t idx) const {
        return bitvector.rank_symbol(idx);
    }


    auto all_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, Sigma>, std::array<uint64_t, Sigma>> {
        auto [rs, prs] = bitvector.all_ranks_and_prefix_ranks(idx);
        return {rs, prs};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bitvector);
    }
};


}
namespace epr8V3 {
template <uint64_t TSigma>
struct OccTable : eprV3_impl::OccTable<TSigma, uint8_t> {
    using eprV3_impl::OccTable<TSigma, uint8_t>::OccTable;
    static auto name() -> std::string {
        return "EPRV3 (8bit)";
    }

    static auto extension() -> std::string {
        return "epr8v3";
    }
};
static_assert(checkOccTable<OccTable>);
}
namespace epr16V3 {
template <uint64_t TSigma>
struct OccTable : eprV3_impl::OccTable<TSigma, uint16_t> {
    using eprV3_impl::OccTable<TSigma, uint16_t>::OccTable;
    static auto name() -> std::string {
        return "EPRV3 (16bit)";
    }

    static auto extension() -> std::string {
        return "epr16v3";
    }
};
static_assert(checkOccTable<OccTable>);
}
namespace epr32V3 {
template <uint64_t TSigma>
struct OccTable : eprV3_impl::OccTable<TSigma, uint32_t> {
    using eprV3_impl::OccTable<TSigma, uint32_t>::OccTable;
    static auto name() -> std::string {
        return "EPRV3 (32bit)";
    }

    static auto extension() -> std::string {
        return "epr32v3";
    }
};
static_assert(checkOccTable<OccTable>);
}

}
}
