#pragma once

#include "concepts.h"

#include <array>
#include <bitset>
#include <cstdint>
#include <vector>


namespace fmindex_collection {
namespace occtable {
namespace interleavedEPRV2_impl {

// counts how many bits are needed to represent the number y
constexpr inline size_t bits_count(size_t y) {
    if (y == 0) return 1;
    size_t i{0};
    while (y != 0) {
        y = y >> 1;
        ++i;
    }
    return i;
}

// computes b to the power of y
constexpr inline size_t pow(size_t b, size_t y) {
    if (y == 0) return 1;
    return pow(b, (y-1)) * b;
}


template <size_t TSigma, size_t TAlignment, typename block_t>
struct Bitvector {

    // number of full length bitvectors needed `2^bitct ≥ TSigma`
    static constexpr auto bitct = bits_count(TSigma-1);
    // next full power of 2
    static constexpr auto bvct  = pow(2, bitct);

    struct alignas(TAlignment) Block {
        std::array<block_t, TSigma> blocks{};
        std::array<uint64_t, bitct> bits{};

        void prefetch() const {
            __builtin_prefetch(reinterpret_cast<void const*>(&blocks), 0, 0);
//            __builtin_prefetch((const void*)(&bits), 0, 0);
        }

        uint64_t rank(size_t idx, size_t symb) const {
            assert(idx >= 0 && idx < 64);
            auto f = [&]<size_t I>(std::index_sequence<I>) {
                return bits[I] ^ -((~symb>>I)&1);
            };
            auto mask = [&]<size_t ...Is>(std::index_sequence<Is...>) {
                return (f(std::index_sequence<Is>{})&...);
            }(std::make_index_sequence<bitct>{});

            auto bitset = std::bitset<64>(mask) << (64-idx);
            return blocks[symb] + bitset.count();
        }

        uint64_t prefix_rank(size_t idx, size_t symb) const {
            auto f = [&]<size_t I>(std::index_sequence<I>, size_t _symb) {
                return bits[I] ^ -((~_symb>>I)&1);
            };
            size_t mask{};

            for (size_t i{0}; i <= symb; ++i) {
                mask |= [&]<size_t ...Is>(std::index_sequence<Is...>) {
                    return (f(std::index_sequence<Is>{}, i)&...);
                }(std::make_index_sequence<bitct>{});
            }

            auto bitset = std::bitset<64>(mask) << (64-idx);
            auto ct = bitset.count();

            for (size_t i{0}; i <= symb; ++i) {
                ct += blocks[i];
            }
            return ct;
        }

        auto all_ranks(size_t idx) const -> std::array<uint64_t, TSigma> {
            assert(idx >= 0 && idx < 64);

            std::array<uint64_t, TSigma> rs{0};

            auto f = [&]<size_t I>(uint64_t symb, std::index_sequence<I>) {
                return bits[I] ^ -((~symb>>I)&1);
            };

            for (size_t i{0}; i < TSigma; ++i) {
                auto mask = [&]<size_t ...Is>(std::index_sequence<Is...>) {
                    return (f(i, std::index_sequence<Is>{})&...);
                }(std::make_index_sequence<bitct>{});
                rs[i] = (std::bitset<64>(mask) << (64 - idx)).count() + blocks[i];
            }
            return rs;
        }

        size_t symbol(size_t idx) const {
            uint64_t symb{};
            for (size_t i{bitct}; i > 0; --i) {
                auto b = (bits[i-1] >> idx) & 1;
                symb = (symb<<1) | b;
            }
            return symb;
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar(blocks, bits);
        }
    };

    static constexpr size_t block_size = sizeof(block_t) * 8;

    std::vector<Block> blocks;
    std::vector<std::array<uint64_t, TSigma>> superBlocks;
    std::array<uint64_t, TSigma+1> C;

    template <typename CB>
    Bitvector(size_t length, CB cb) {
        blocks.reserve(length/64+2);

        std::array<uint64_t, TSigma> sblock_acc{0};
        std::array<block_t, TSigma> block_acc{0};

        for (size_t size{0}; size < length; ++size) {
            if (size % (1ul<<block_size) == 0) { // new super block + new block
                superBlocks.emplace_back(sblock_acc);
                blocks.emplace_back();
                block_acc = {};
            } else if (size % 64 == 0) { // new block
                blocks.emplace_back();
                blocks.back().blocks = block_acc;
            }
            auto blockId      = size >>  6;
            auto bitId        = size &  63;

            size_t symb = cb(size);

            for (size_t i{}; i < bitct; ++i) {
                auto b = ((symb>>i)&1);
                blocks[blockId].bits[i] |= (b << bitId);
            }
            block_acc[symb] += 1;
            sblock_acc[symb] += 1;
        }

        C[0] = 0;
        for (size_t i{0}; i < TSigma; ++i) {
            C[i+1] = sblock_acc[i] + C[i];
        }
    }

    Bitvector(cereal_tag) {}


    size_t memoryUsage() const {
        return blocks.size() * sizeof(blocks.back())
            + superBlocks.size() * sizeof(superBlocks.back())
            + sizeof(C);
    }

    void prefetch(size_t idx) const {
        auto blockId      = idx >>  6;
        blocks[blockId].prefetch();
    }

    uint64_t rank(uint64_t idx, size_t symb) const {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;
        auto bitId        = idx &  63;
        return blocks[blockId].rank(bitId, symb) + superBlocks[superBlockId][symb] + C[symb];
    }

    uint64_t prefix_rank(uint64_t idx, size_t symb) const {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;
        auto bitId        = idx &  63;
        uint64_t a={};
        for (size_t i{0}; i<= symb; ++i) {
            a += superBlocks[superBlockId][i];
        }
        return blocks[blockId].prefix_rank(bitId, symb) + a;
    }


    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> block_size;
        auto bitId        = idx &  63;
        auto res = std::array<uint64_t, TSigma>{};
        for (size_t symb{0}; symb < TSigma; ++symb) {
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
        for (size_t symb{1}; symb < TSigma; ++symb) {
            prs[symb] = prs[symb-1] + superBlocks[superBlockId][symb] + rs[symb];
            rs[symb] += C[symb] + superBlocks[superBlockId][symb];
        }
        return {rs, prs};
    }

    size_t symbol(uint64_t idx) const {
        auto blockId      = idx >>  6;
        auto bitId        = idx &  63;
        return blocks[blockId].symbol(bitId);
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(blocks, superBlocks, C);
    }
};


template <size_t TSigma, typename block_t, size_t TAlignment>
struct OccTable {
    using TLengthType = uint64_t;
    static constexpr size_t Sigma = TSigma;

    Bitvector<Sigma, TAlignment, block_t> bitvector;

    static size_t expectedMemoryUsage(size_t length) {
        using Block = typename Bitvector<TSigma, TAlignment, block_t>::Block;
        auto blockSize = std::max(alignof(Block), sizeof(Block));

        size_t C           = sizeof(uint64_t) * (Sigma+1);
        size_t blocks      = blockSize        * (length+1) / 64;
        size_t superblocks = sizeof(uint64_t) * (length+1) / (1ul << (sizeof(block_t) * 8));
        return C + blocks + superblocks;
    }

    OccTable(std::vector<uint8_t> const& _bwt)
        : bitvector(_bwt.size(), [&](size_t i) -> uint8_t {
            return _bwt[i];
        })
    {}

    OccTable(cereal_tag)
        : bitvector(cereal_tag{})
    {}

    size_t memoryUsage() const {
        return bitvector.memoryUsage() + sizeof(OccTable);
    }

    uint64_t size() const {
        return bitvector.C.back();
    }

    auto prefetch(uint64_t idx) const {
        bitvector.prefetch(idx);
    }

    uint64_t rank(uint64_t idx, size_t symb) const {
        return bitvector.rank(idx, symb);
    }

    uint64_t prefix_rank(uint64_t idx, size_t symb) const {
        return bitvector.prefix_rank(idx, symb);
    }

    size_t symbol(uint64_t idx) const {
        return bitvector.symbol(idx);
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
template <size_t TSigma>
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
template <size_t TSigma>
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
template <size_t TSigma>
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
template <size_t TSigma>
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
template <size_t TSigma>
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
template <size_t TSigma>
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
