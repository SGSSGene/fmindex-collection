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
#if __has_include(<cereal/archives/binary.hpp>)
#include <cereal/archives/binary.hpp>
#endif
#include <cstdint>
#include <cstring>
#include <span>
#include <stdexcept>
#include <vector>


namespace fmindex_collection {
namespace occtable {
namespace interleavedEPRV2_impl {


template <size_t I, size_t N>
auto requested_bit_for_symbol(uint64_t symb, std::array<uint64_t, N> const& bits) -> uint64_t {
    // I bit of symb is symbol of interest
    // 1. Detect inversed I bit
    auto inversed_bit = (~symb>>I) & 1;
    // 2. Expand 0 → 000...000    and 1 → 111...111
    auto expanded_bit = - inversed_bit;
    // 3. Return all bits that are not set as set in expanded_bit
    return bits[I] ^ expanded_bit;
}


template <size_t N>
auto bits_have_symbols(uint64_t symb, std::array<uint64_t, N> const& bits) -> uint64_t {
    auto bits_set = [&]<uint64_t ...Is>(std::integer_sequence<uint64_t, Is...>) {
        return (requested_bit_for_symbol<Is>(symb, bits)&...);
    }(std::make_integer_sequence<uint64_t, N>{});
    return bits_set;
}

template <size_t N>
auto bits_have_symbols_or_less(uint64_t symb, std::array<uint64_t, N> const& bits) -> uint64_t {
    size_t bit_set{};
    for (uint64_t i{0}; i <= symb; ++i) {
        bit_set |= bits_have_symbols(i, bits);
    }
    return bit_set;
}


template <size_t Sigma, size_t N>
auto bits_have_symbols_or_less2(uint64_t symb, std::array<uint64_t, N> const& bits) -> uint64_t {
    size_t bit_set{};
    for (uint64_t i{0}; i < Sigma; ++i) {
        bit_set |= bits_have_symbols(symb, bits) * (size_t)(i <= symb);
    }
    return bit_set;
}


template <uint64_t TSigma, uint64_t TAlignment, typename block_t>
struct Bitvector {

    // number of full length bitvectors needed `2^sigma_bits ≥ TSigma`
    static constexpr auto sigma_bits = required_bits(TSigma-1);

    struct alignas(TAlignment) Block {
        std::array<block_t, TSigma> blocks{};
        std::array<uint64_t, sigma_bits> bits{};

        void prefetch() const {
            __builtin_prefetch(reinterpret_cast<void const*>(&blocks), 0, 0);
//            __builtin_prefetch((const void*)(&bits), 0, 0);
        }

        uint64_t rank(uint64_t idx, uint64_t symb) const {
            assert(idx < 64);
            auto bits_set = bits_have_symbols(symb, bits);
            auto bits_masked = std::bitset<64>(bits_set) << (64-idx);
            return blocks[symb] + bits_masked.count();
        }

        uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
            auto bits_set = bits_have_symbols_or_less(symb, bits);

            auto bits_masked = std::bitset<64>(bits_set) << (64-idx);
            auto ct = bits_masked.count();

            for (uint64_t i{0}; i <= symb; ++i) {
                ct += blocks[i];
            }
            return ct;
        }

        auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
            assert(idx < 64);

            std::array<uint64_t, TSigma> rs;

            for (uint64_t i{0}; i < TSigma; ++i) {
                auto bits_set = bits_have_symbols(i, bits);
                rs[i] = (std::bitset<64>(bits_set) << (64 - idx)).count() + blocks[i];
            }
            return rs;
        }

        uint64_t symbol(uint64_t idx) const {
            uint64_t symb{};
            for (uint64_t i{sigma_bits}; i > 0; --i) {
                auto b = (bits[i-1] >> idx) & 1;
                symb = (symb<<1) | b;
            }
            return symb;
        }

        auto rank_symbol(uint64_t idx) const -> std::tuple<uint64_t, uint64_t> {
            assert(idx < 64);

            uint64_t symb{};
            uint64_t mask{};
            auto f = [&]<auto I>(std::integer_sequence<uint64_t, I>) {
                auto b = (bits[uint64_t{I}] >> idx) & 1;
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

    static constexpr uint64_t block_size = sizeof(block_t) * 8;

    std::vector<Block> blocks;
    std::vector<std::array<uint64_t, TSigma>> superBlocks;
    std::array<uint64_t, TSigma+1> C;

    Bitvector(std::span<uint8_t const> _bwt) {
        // Next three lines are a reserve call, with zero initialization
        // This is required, so padding bytes will also be zero
        blocks.resize(_bwt.size()/64+1);
        memset((void*)blocks.data(), 0, blocks.size() * sizeof(Block));
        blocks.resize(0);

        auto sblock_acc = std::array<uint64_t, TSigma>{}; // accumulator for super blocks
        auto block_acc  = std::array<block_t, TSigma>{};  // accumulator for blocks

        for (uint64_t size{0}; size < _bwt.size();) {
            superBlocks.emplace_back(sblock_acc);
            block_acc = {};

            for (uint64_t blockId{0}; blockId < (1ull<<block_size)/64 and size < _bwt.size(); ++blockId) {
                blocks.emplace_back();
                blocks.back().blocks = block_acc;

                for (uint64_t bitId{0}; bitId < 64 and size < _bwt.size(); ++bitId, ++size) {

                    uint64_t symb = _bwt[size];

                    for (uint64_t i{}; i < sigma_bits; ++i) {
                        auto b = ((symb>>i)&1);
                        blocks.back().bits[i] |= (b << bitId);
                    }

                    block_acc[symb] += 1;
                    sblock_acc[symb] += 1;
                }
            }
        }
        // Add a new block, so we can access one row more than our bwt length
        if (_bwt.size() % 64 == 0) {
            superBlocks.emplace_back(sblock_acc);
            blocks.emplace_back();
            blocks.back().blocks = block_acc;
        }

        C[0] = 0;
        for (uint64_t i{0}; i < TSigma; ++i) {
            C[i+1] = sblock_acc[i] + C[i];
        }
    }

    Bitvector(cereal_tag) {}


    uint64_t memoryUsage() const {
        return blocks.size() * sizeof(blocks.back())
            + superBlocks.size() * sizeof(superBlocks.back())
            + sizeof(C);
    }

    void prefetch(uint64_t idx) const {
        auto blockId      = idx >>  6;
        blocks[blockId].prefetch();
    }

    uint64_t rank(uint64_t idx, uint64_t symb) const {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;
        auto bitId        = idx &  63;
        return blocks[blockId].rank(bitId, symb) + superBlocks[superBlockId][symb] + C[symb];
    }

    uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;
        auto bitId        = idx &  63;
        uint64_t a={};
        for (uint64_t i{0}; i<= symb; ++i) {
            a += superBlocks[superBlockId][i];
        }
        return blocks[blockId].prefix_rank(bitId, symb) + a;
    }


    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;
        auto bitId        = idx &  63;
        auto res = std::array<uint64_t, TSigma>{};
        for (uint64_t symb{0}; symb < TSigma; ++symb) {
            res[symb] = blocks[blockId].rank(bitId, symb) + superBlocks[superBlockId][symb] + C[symb];
        }
        return res;
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, TSigma>, std::array<uint64_t, TSigma>> {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;
        auto bitId        = idx &  63;

        std::array<uint64_t, TSigma> prs;
        auto rs = blocks[blockId].all_ranks(bitId);

        rs[0] += superBlocks[superBlockId][0] + C[0];
        prs[0] = rs[0];
        for (uint64_t symb{1}; symb < TSigma; ++symb) {
            prs[symb] = prs[symb-1] + superBlocks[superBlockId][symb] + rs[symb];
            rs[symb] += C[symb] + superBlocks[superBlockId][symb];
        }
        return {rs, prs};
    }

    uint64_t symbol(uint64_t idx) const {
        auto blockId      = idx >>  6;
        auto bitId        = idx &  63;
        return blocks[blockId].symbol(bitId);
    }

    uint64_t rank_symbol(uint64_t idx) const {
        auto blockId = idx >> 6;
        auto bitId = idx & 63;
        auto superBlockId = idx >> block_size;
        auto [rank, symb] = blocks[blockId].rank_symbol(bitId);
        return rank + superBlocks[superBlockId][symb] + C[symb];
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
            ar(blocks, superBlocks, C);
        } else if (version == 1) {
#if __has_include(<cereal/archives/binary.hpp>)
            if constexpr (std::same_as<Archive, cereal::BinaryOutputArchive>
                            || std::same_as<Archive, cereal::BinaryInputArchive>) {
                auto l = blocks.size();
                ar(l);
                blocks.resize(l);
                ar(cereal::binary_data(blocks.data(), l * sizeof(Block)),
                   superBlocks,
                   C);
            } else
#endif
            throw std::runtime_error("fmindex-collection - InterleavedEPRV2 was created with binary data, but this is not available in this app");
        } else {
            throw std::runtime_error("fmindex-collection - InterleavedEPRV2 was created with legacy format - not readable by this app");
        }

    }
};


template <uint64_t TSigma, typename block_t, uint64_t TAlignment>
struct OccTable {
    using TLengthType = uint64_t;
    static constexpr uint64_t Sigma = TSigma;

    Bitvector<Sigma, TAlignment, block_t> bitvector;

    static uint64_t expectedMemoryUsage(uint64_t length) {
        using Block = typename Bitvector<TSigma, TAlignment, block_t>::Block;
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
namespace interleavedEPR8V2 {
template <uint64_t TSigma>
struct OccTable : interleavedEPRV2_impl::OccTable<TSigma, uint8_t, 8> {
    using interleavedEPRV2_impl::OccTable<TSigma, uint8_t, 8>::OccTable;
    static auto name() -> std::string {
        return "Interleaved EPRV2 (8bit)";
    }

    static auto extension() -> std::string {
        return "iepr8v2";
    }
};
static_assert(checkOccTable<OccTable>);
}
namespace interleavedEPR16V2 {
template <uint64_t TSigma>
struct OccTable : interleavedEPRV2_impl::OccTable<TSigma, uint16_t, 8> {
    using interleavedEPRV2_impl::OccTable<TSigma, uint16_t, 8>::OccTable;
    static auto name() -> std::string {
        return "Interleaved EPRV2 (16bit)";
    }

    static auto extension() -> std::string {
        return "iepr16v2";
    }
};
static_assert(checkOccTable<OccTable>);
}
namespace interleavedEPR32V2 {
template <uint64_t TSigma>
struct OccTable : interleavedEPRV2_impl::OccTable<TSigma, uint32_t, 8> {
    using interleavedEPRV2_impl::OccTable<TSigma, uint32_t, 8>::OccTable;
    static auto name() -> std::string {
        return "Interleaved EPRV2 (32bit)";
    }

    static auto extension() -> std::string {
        return "iepr32v2";
    }
};
static_assert(checkOccTable<OccTable>);
}
namespace interleavedEPR8V2Aligned {
template <uint64_t TSigma>
struct OccTable : interleavedEPRV2_impl::OccTable<TSigma, uint8_t, 64> {
    using interleavedEPRV2_impl::OccTable<TSigma, uint8_t, 64>::OccTable;
    static auto name() -> std::string {
        return "Interleaved EPRV2 (8bit, aligned)";
    }

    static auto extension() -> std::string {
        return "iepr8v2a";
    }
};
static_assert(checkOccTable<OccTable>);
}
namespace interleavedEPR16V2Aligned {
template <uint64_t TSigma>
struct OccTable : interleavedEPRV2_impl::OccTable<TSigma, uint16_t, 64> {
    using interleavedEPRV2_impl::OccTable<TSigma, uint16_t, 64>::OccTable;
    static auto name() -> std::string {
        return "Interleaved EPRV2 (16bit, aligned)";
    }

    static auto extension() -> std::string {
        return "iepr16v2a";
    }
};
static_assert(checkOccTable<OccTable>);
}
namespace interleavedEPR32V2Aligned {
template <uint64_t TSigma>
struct OccTable : interleavedEPRV2_impl::OccTable<TSigma, uint32_t, 64> {
    using interleavedEPRV2_impl::OccTable<TSigma, uint32_t, 64>::OccTable;
    static auto name() -> std::string {
        return "Interleaved EPRV2 (32bit, aligned)";
    }

    static auto extension() -> std::string {
        return "iepr32v2a";
    }
};
static_assert(checkOccTable<OccTable>);

}

}
}
