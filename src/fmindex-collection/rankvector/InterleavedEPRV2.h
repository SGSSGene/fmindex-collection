// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../builtins.h"
#include "concepts.h"
#include "utils.h"

#include <bitset>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <span>
#include <stdexcept>
#include <vector>

#if __has_include(<cereal/archives/binary.hpp>)
#include <cereal/archives/binary.hpp>
#endif

namespace fmindex_collection::rankvector {

template <size_t TSigma, uint64_t TAlignment, typename block_t>
struct InterleavedEPRV2 {
    static constexpr size_t Sigma = TSigma;

    template <size_t I, size_t N>
    static auto requested_bit_for_symbol(uint64_t symb, std::array<uint64_t, N> const& bits) -> uint64_t {
        // I bit of symb is symbol of interest
        // 1. Detect inversed I bit
        auto inversed_bit = (~symb>>I) & 1;
        // 2. Expand 0 → 000...000    and 1 → 111...111
        auto expanded_bit = - inversed_bit;
        // 3. Return all bits that are not set as set in expanded_bit
        return bits[I] ^ expanded_bit;
    }


    template <size_t N>
    static auto bits_have_symbols(uint64_t symb, std::array<uint64_t, N> const& bits) -> uint64_t {
        auto bits_set = [&]<uint64_t ...Is>(std::integer_sequence<uint64_t, Is...>) {
            return (requested_bit_for_symbol<Is>(symb, bits)&...);
        }(std::make_integer_sequence<uint64_t, N>{});
        return bits_set;
    }

    template <size_t N>
    static auto bits_have_symbols_or_less(uint64_t symb, std::array<uint64_t, N> const& bits) -> uint64_t {
        uint64_t bit_set{};
        for (size_t i{0}; i <= symb; ++i) {
            bit_set |= bits_have_symbols(i, bits);
        }
        return bit_set;
    }

    // number of full length bitvectors needed `2^sigma_bits ≥ TSigma`
    static constexpr auto sigma_bits = required_bits(TSigma-1);

    struct alignas(TAlignment) Block {
        std::array<block_t, TSigma> blocks{};
        std::array<uint64_t, sigma_bits> bits{};

        void prefetch() const {
            __builtin_prefetch(reinterpret_cast<void const*>(&blocks), 0, 0);
//            __builtin_prefetch((const void*)(&bits), 0, 0);
        }

        uint64_t rank(size_t idx, uint64_t symb) const {
            assert(idx < 64);
            auto bits_set = bits_have_symbols(symb, bits);
            auto bits_masked = std::bitset<64>(bits_set) << (64-idx);
            return blocks[symb] + bits_masked.count();
        }

        uint64_t prefix_rank(size_t idx, uint64_t symb) const {
            auto bits_set = bits_have_symbols_or_less(symb, bits);

            auto bits_masked = std::bitset<64>(bits_set) << (64-idx);
            auto ct = bits_masked.count();

            for (size_t i{0}; i <= symb; ++i) {
                ct += blocks[i];
            }
            return ct;
        }

        auto all_ranks(size_t idx) const -> std::array<uint64_t, TSigma> {
            assert(idx < 64);

            auto rs = std::array<uint64_t, TSigma>{};

            for (size_t i{0}; i < TSigma; ++i) {
                auto bits_set = bits_have_symbols(i, bits);
                rs[i] = (std::bitset<64>(bits_set) << (64 - idx)).count() + blocks[i];
            }
            return rs;
        }

        uint64_t symbol(size_t idx) const {
            uint64_t symb{};
            for (size_t i{sigma_bits}; i > 0; --i) {
                auto b = (bits[i-1] >> idx) & 1;
                symb = (symb<<1) | b;
            }
            return symb;
        }

        auto rank_symbol(size_t idx) const -> std::tuple<uint64_t, uint64_t> {
            assert(idx < 64);

            uint64_t symb{};
            uint64_t mask{};
            auto f = [&]<uint64_t I>(std::integer_sequence<uint64_t, I>) {
                auto b = (bits[I] >> idx) & 1;
                mask |= bits[I] ^ -b;
                symb |= b << I;
            };
            [&]<auto ...Is>(std::integer_sequence<uint64_t, Is...>) {
                (f(std::integer_sequence<uint64_t, Is>{}) ,...);
            }(std::make_integer_sequence<uint64_t, sigma_bits>{});

            auto bitset = std::bitset<64>{~mask} << (64-idx);

            return {blocks[symb] + bitset.count(), symb};
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar(blocks, bits);
        }
    };

    static constexpr size_t block_size = sizeof(block_t) * 8;

    std::vector<Block> blocks;
    std::vector<std::array<uint64_t, TSigma>> superBlocks;
    size_t totalLength{};

    InterleavedEPRV2() = default;
    InterleavedEPRV2(std::span<uint8_t const> _symbols) {
        totalLength = _symbols.size();
        // Next three lines are a reserve call, with zero initialization
        // This is required, so padding bytes will also be zero
        blocks.resize(_symbols.size()/64+1);
        memset((void*)blocks.data(), 0, blocks.size() * sizeof(Block));
        blocks.resize(0);

        auto sblock_acc = std::array<uint64_t, TSigma>{}; // accumulator for super blocks
        auto block_acc  = std::array<block_t, TSigma>{};  // accumulator for blocks

        for (size_t size{0}; size < _symbols.size();) {
            superBlocks.emplace_back(sblock_acc);
            block_acc = {};

            for (size_t blockId{0}; blockId < (uint64_t{1}<<block_size)/64 and size < _symbols.size(); ++blockId) {
                blocks.emplace_back();
                blocks.back().blocks = block_acc;

                for (size_t bitId{0}; bitId < 64 and size < _symbols.size(); ++bitId, ++size) {

                    uint64_t symb = _symbols[size];

                    for (size_t i{}; i < sigma_bits; ++i) {
                        auto b = ((symb>>i)&1);
                        blocks.back().bits[i] |= (b << bitId);
                    }

                    block_acc[symb] += 1;
                    sblock_acc[symb] += 1;
                }
            }
        }
        // Add a new block, so we can access one row more than our symbols length
        if (_symbols.size() % 64 == 0) {
            superBlocks.emplace_back(sblock_acc);
            blocks.emplace_back();
            blocks.back().blocks = block_acc;
        }
    }

    size_t size() const {
        return totalLength;
    }

    void prefetch(size_t idx) const {
        auto blockId      = idx >>  6;
        blocks[blockId].prefetch();
    }

    uint8_t symbol(size_t idx) const {
        auto blockId      = idx >>  6;
        auto bitId        = idx &  63;
        return blocks[blockId].symbol(bitId);
    }

    uint64_t rank(size_t idx, uint8_t symb) const {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;
        auto bitId        = idx &  63;
        return blocks[blockId].rank(bitId, symb) + superBlocks[superBlockId][symb];
    }

    uint64_t prefix_rank(size_t idx, uint8_t symb) const {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;
        auto bitId        = idx &  63;
        uint64_t a{};
        for (size_t i{0}; i<= symb; ++i) {
            a += superBlocks[superBlockId][i];
        }
        return blocks[blockId].prefix_rank(bitId, symb) + a;
    }


    auto all_ranks(size_t idx) const -> std::array<uint64_t, TSigma> {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;
        auto bitId        = idx &  63;
        auto res = std::array<uint64_t, TSigma>{};
        for (size_t symb{0}; symb < TSigma; ++symb) {
            res[symb] = blocks[blockId].rank(bitId, symb) + superBlocks[superBlockId][symb];
        }
        return res;
    }

    auto all_ranks_and_prefix_ranks(size_t idx) const -> std::tuple<std::array<uint64_t, TSigma>, std::array<uint64_t, TSigma>> {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;
        auto bitId        = idx &  63;

        auto prs = std::array<uint64_t, TSigma>{};
        auto rs = blocks[blockId].all_ranks(bitId);

        rs[0] += superBlocks[superBlockId][0];
        prs[0] = rs[0];
        for (size_t symb{1}; symb < TSigma; ++symb) {
            prs[symb] = prs[symb-1] + superBlocks[superBlockId][symb] + rs[symb];
            rs[symb] += superBlocks[superBlockId][symb];
        }
        return {rs, prs};
    }

    uint64_t rank_symbol(size_t idx) const {
        auto blockId      = idx >> 6;
        auto bitId        = idx & 63;
        auto superBlockId = idx >> block_size;
        auto [rank, symb] = blocks[blockId].rank_symbol(bitId);
        return rank + superBlocks[superBlockId][symb];
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        // 0 version: slow path
        // 1 version: binary path (fast)
        auto version = []() {
#if __has_include(<cereal/archives/binary.hpp>)
            if constexpr (std::same_as<Archive, cereal::BinaryOutputArchive>
                            || std::same_as<Archive, cereal::BinaryInputArchive>) {
                return 1;
            }
#endif
            return 0;
        }();
        ar(version);

        if (version == 0) {
            ar(blocks, superBlocks, totalLength);
        } else if (version == 1) {
#if __has_include(<cereal/archives/binary.hpp>)
            if constexpr (std::same_as<Archive, cereal::BinaryOutputArchive>
                            || std::same_as<Archive, cereal::BinaryInputArchive>) {
                auto l = blocks.size();
                ar(l);
                blocks.resize(l);
                ar(cereal::binary_data(blocks.data(), l * sizeof(Block)),
                   superBlocks,
                   totalLength);
            } else
#endif
            throw std::runtime_error("fmindex-collection - InterleavedEPRV2 was created with binary data, but this is not available in this app");
        } else {
            throw std::runtime_error("fmindex-collection - InterleavedEPRV2 was created with legacy format - not readable by this app");
        }

    }

};

template <size_t TSigma> using InterleavedEPRV2_8         = InterleavedEPRV2<TSigma,  8,  uint8_t>;
static_assert(checkRankVector<InterleavedEPRV2_8>);

template <size_t TSigma> using InterleavedEPRV2_16        = InterleavedEPRV2<TSigma,  8, uint16_t>;
static_assert(checkRankVector<InterleavedEPRV2_16>);

#if SIZE_MAX == UINT64_MAX
template <size_t TSigma> using InterleavedEPRV2_32        = InterleavedEPRV2<TSigma,  8, uint32_t>;
static_assert(checkRankVector<InterleavedEPRV2_32>);
#endif

template <size_t TSigma> using InterleavedEPRV2_8Aligned  = InterleavedEPRV2<TSigma, 64,  uint8_t>;
static_assert(checkRankVector<InterleavedEPRV2_8Aligned>);

template <size_t TSigma> using InterleavedEPRV2_16Aligned = InterleavedEPRV2<TSigma, 64, uint16_t>;
static_assert(checkRankVector<InterleavedEPRV2_16Aligned>);

#if SIZE_MAX == UINT64_MAX
template <size_t TSigma> using InterleavedEPRV2_32Aligned = InterleavedEPRV2<TSigma, 64, uint32_t>;
static_assert(checkRankVector<InterleavedEPRV2_32Aligned>);
#endif



}
