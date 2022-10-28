#pragma once

#include "concepts.h"

#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#if __has_include(<cereal/archives/binary.hpp>)
#include <cereal/archives/binary.hpp>
#endif
#include <cstdint>
#include <vector>


namespace fmindex_collection {
namespace occtable {
namespace interleavedEPRV2_impl {

// counts how many bits are needed to represent the number y
constexpr inline uint64_t bits_count(uint64_t y) {
    if (y == 0) return 1;
    uint64_t i{0};
    while (y != 0) {
        y = y >> 1;
        ++i;
    }
    return i;
}

// computes b to the power of y
constexpr inline uint64_t pow(uint64_t b, uint64_t y) {
    if (y == 0) return 1;
    return pow(b, (y-1)) * b;
}


template <uint64_t TSigma, uint64_t TAlignment, typename block_t>
struct Bitvector {

    // number of full length bitvectors needed `2^bitct ≥ TSigma`
    static constexpr auto bitct = bits_count(TSigma-1);
    // next full power of 2
    static constexpr auto bvct  = pow(2, bitct);

    struct alignas(TAlignment) Block {
        std::array<block_t, TSigma> blocks;
        std::array<uint64_t, bitct> bits;

        void prefetch() const {
            __builtin_prefetch(reinterpret_cast<void const*>(&blocks), 0, 0);
//            __builtin_prefetch((const void*)(&bits), 0, 0);
        }

        uint64_t rank(uint64_t idx, uint64_t symb) const {
            assert(idx < 64);
            auto f = [&]<uint64_t I>(std::integer_sequence<uint64_t, I>) {
                return bits[I] ^ -((~symb>>I)&1);
            };
            auto mask = [&]<uint64_t ...Is>(std::integer_sequence<uint64_t, Is...>) {
                return (f(std::integer_sequence<uint64_t, Is>{})&...);
            }(std::make_integer_sequence<uint64_t, bitct>{});

            auto bitset = std::bitset<64>(mask) << (64-idx);
            return blocks[symb] + bitset.count();
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
            auto ct = bitset.count();

            for (uint64_t i{0}; i <= symb; ++i) {
                ct += blocks[i];
            }
            return ct;
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
                rs[i] = (std::bitset<64>(mask) << (64 - idx)).count() + blocks[i];
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
            auto f = [&]<auto I>(std::integer_sequence<uint64_t, I>) {
                auto b = (bits[uint64_t{I}] >> idx) & 1;
                mask |= bits[I] ^ -b;
                symb |= b << I;
            };
            [&]<auto ...Is>(std::integer_sequence<uint64_t, Is...>) {
                (f(std::integer_sequence<uint64_t, Is>{}) ,...);
            }(std::make_integer_sequence<uint64_t, bitct>{});

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

    template <typename CB>
    Bitvector(uint64_t length, CB cb) {
        blocks.reserve(length/64+2);

        std::array<uint64_t, TSigma> sblock_acc{0};
        std::array<block_t, TSigma> block_acc{0};

        for (uint64_t size{0}; size < length; ++size) {
            if (size % (1ul<<block_size) == 0) { // new super block + new block
                superBlocks.emplace_back(sblock_acc);
                blocks.emplace_back(Block{});
                block_acc = {};
            } else if (size % 64 == 0) { // new block
                blocks.emplace_back(Block{});
                blocks.back().blocks = block_acc;
            }
            auto blockId      = size >>  6;
            auto bitId        = size &  63;

            uint64_t symb = cb(size);

            for (uint64_t i{}; i < bitct; ++i) {
                auto b = ((symb>>i)&1);
                blocks[blockId].bits[i] |= (b << bitId);
            }
            block_acc[symb] += 1;
            sblock_acc[symb] += 1;
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
        {
            ar(blocks, superBlocks, C);
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
        uint64_t superblocks = sizeof(uint64_t) * (length+1) / (1ul << (sizeof(block_t) * 8));
        return C + blocks + superblocks;
    }

    OccTable(std::vector<uint8_t> const& _bwt)
        : bitvector(_bwt.size(), [&](uint64_t i) -> uint8_t {
            return _bwt[i];
        })
    {}

    OccTable(cereal_tag)
        : bitvector(cereal_tag{})
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
struct OccTable : interleavedEPRV2_impl::OccTable<TSigma, uint8_t, 1> {
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
struct OccTable : interleavedEPRV2_impl::OccTable<TSigma, uint16_t, 1> {
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
struct OccTable : interleavedEPRV2_impl::OccTable<TSigma, uint32_t, 1> {
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
