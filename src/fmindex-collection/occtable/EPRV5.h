#pragma once

#include "concepts.h"

#include <array>
#include <bitset>
#include <cstdint>
#include <vector>


namespace fmindex_collection {
namespace occtable {
namespace eprV5_impl {

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


template <size_t TSigma, size_t TAlignment>
struct Bitvector {

    // number of full length bitvectors needed `2^bitct ≥ TSigma`
    static constexpr auto bitct = bits_count(TSigma-1);
    // next full power of 2
    static constexpr auto bvct  = pow(2, bitct);

    struct InBits {
        std::array<uint64_t, bitct> bits{};

        uint64_t rank(size_t idx, size_t symb) const {
            assert(idx < 64);
            auto f = [&]<size_t I>(std::index_sequence<I>) {
                return bits[I] ^ -((~symb>>I)&1);
            };
            auto mask = [&]<size_t ...Is>(std::index_sequence<Is...>) {
                return (f(std::index_sequence<Is>{})&...);
            }(std::make_index_sequence<bitct>{});

            auto bitset = std::bitset<64>(mask) << (64-idx);
            return bitset.count();
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
            return bitset.count();
        }

        auto all_ranks(size_t idx) const -> std::array<uint64_t, TSigma> {
            assert(idx < 64);

            std::array<uint64_t, TSigma> rs{0};

            auto f = [&]<size_t I>(uint64_t symb, std::index_sequence<I>) {
                return bits[I] ^ -((~symb>>I)&1);
            };

            for (size_t i{0}; i < TSigma; ++i) {
                auto mask = [&]<size_t ...Is>(std::index_sequence<Is...>) {
                    return (f(i, std::index_sequence<Is>{})&...);
                }(std::make_index_sequence<bitct>{});
                rs[i] = (std::bitset<64>(mask) << (64 - idx)).count();
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

        auto rank_symbol(size_t idx) const -> std::tuple<size_t, size_t> {
            assert(idx < 64);

            uint64_t symb{};
            uint64_t mask{};
            auto f = [&]<size_t I>(std::index_sequence<I>) {
                auto b = (bits[I] >> idx) & 1;
                mask |= bits[I] ^ -b;
                symb |= b << I;
            };
            [&]<size_t ...Is>(std::index_sequence<Is...>) {
                (f(std::index_sequence<Is>{}) ,...);
            }(std::make_index_sequence<bitct>{});

            auto bitset = std::bitset<64>{~mask} << (64-idx);

            return {bitset.count(), symb};
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar(bits);
        }
    };

    using blockL0_t = uint8_t;
    using blockL1_t = uint16_t;
//    using blockL2_t = uint32_t;
    static constexpr size_t level0_size = sizeof(blockL0_t) * 8;
    static constexpr size_t level1_size = sizeof(blockL1_t) * 8;
//    static constexpr size_t level2_size = sizeof(blockL2_t) * 8;

    using BlockL0 = std::array<blockL0_t, TSigma>;
    using BlockL1 = std::array<blockL1_t, TSigma>;
//    using BlockL2 = std::array<blockL2_t, TSigma>;

    std::vector<InBits> bits;
    std::vector<BlockL0> level0;
    std::vector<BlockL1> level1;
//    std::vector<BlockL2> level2;
    std::vector<std::array<uint64_t, TSigma>> superBlocks;


    std::array<uint64_t, TSigma+1> C;

    template <typename CB>
    Bitvector(size_t length, CB cb) {
//        level2.reserve(length/(1ul<<level2_size)+2);
        level1.reserve(length/(1ul<<level1_size)+2);
        level0.reserve(length/64+2);
        bits.reserve(length/64+2);

        std::array<blockL0_t, TSigma> blockL0_acc{0};
        std::array<blockL1_t, TSigma> blockL1_acc{0};
//        std::array<blockL2_t, TSigma> blockL2_acc{0};
        std::array<uint64_t, TSigma> sblock_acc{0};


        for (size_t size{0}; size < length; ++size) {
            if (size % (1ul<<level1_size) == 0) { // new l3 block
                superBlocks.emplace_back(sblock_acc);
                level1.emplace_back();
                level0.emplace_back();
                bits.emplace_back();
                blockL0_acc = {};
                blockL1_acc = {};
/*            } else if (size % (1ul<<level1_size) == 0) { // new l2 block
                level2.emplace_back(blockL2_acc);
                level1.emplace_back();
                level0.emplace_back();
                bits.emplace_back();
                blockL0_acc = {};
                blockL1_acc = {};*/
            } else if (size % (1ul<<level0_size) == 0) { // new l1 block
                level1.emplace_back(blockL1_acc);
                level0.emplace_back();
                bits.emplace_back();

                blockL0_acc = {};
            } else if (size % 64 == 0) { // new l0 block
                level0.emplace_back();
                bits.emplace_back();
                level0.back() = blockL0_acc;
            }
            auto level0Id     = size >>  6;
            auto bitId        = size &  63;

            size_t symb = cb(size);

            for (size_t i{}; i < bitct; ++i) {
                auto b = ((symb>>i)&1);
                bits[level0Id].bits[i] |= (b << bitId);
            }
            blockL0_acc[symb] += 1;
            blockL1_acc[symb] += 1;
            sblock_acc[symb] += 1;
        }

        C[0] = 0;
        for (size_t i{0}; i < TSigma; ++i) {
            C[i+1] = sblock_acc[i] + C[i];
        }
    }

    Bitvector(cereal_tag) {}


    size_t memoryUsage() const {
        return    bits.size() * sizeof(bits.back())
                + level0.size() * sizeof(level0.back())
                + level1.size() * sizeof(level1.back())
//                + level2.size() * sizeof(level2.back())
                + superBlocks.size() * sizeof(superBlocks.back())
                + sizeof(C);
    }

    void prefetch(size_t idx) const {
        auto level0Id     = idx >>  6;
        auto level1Id     = idx >> level0_size;
//        auto level2Id     = idx >> level1_size;
        auto superBlockId = idx >> level1_size;

//        __builtin_prefetch(reinterpret_cast<void const*>(&level2[level2Id]), 0, 0);
        __builtin_prefetch(reinterpret_cast<void const*>(&level1[level1Id]), 0, 0);
        __builtin_prefetch(reinterpret_cast<void const*>(&level0[level0Id]), 0, 0);
        __builtin_prefetch(reinterpret_cast<void const*>(&bits[level1Id]), 0, 0);
        __builtin_prefetch(reinterpret_cast<void const*>(&superBlocks[superBlockId]), 0, 0);
    }

    uint64_t rank(uint64_t idx, size_t symb) const {
        prefetch(idx);

        auto level0Id     = idx >>  6;
        auto level1Id     = idx >> level0_size;
//        auto level2Id     = idx >> level1_size;
        auto superBlockId = idx >> level1_size;
        auto bitId        = idx &  63;
        return    bits[level0Id].rank(bitId, symb)
                + level0[level0Id][symb]
                + level1[level1Id][symb]
//                + level2[level2Id][symb]
                + superBlocks[superBlockId][symb]
                + C[symb];
    }

    uint64_t prefix_rank(uint64_t idx, size_t symb) const {
        prefetch(idx);

        auto level0Id     = idx >>  6;
        auto level1Id     = idx >> level0_size;
//        auto level2Id     = idx >> level1_size;
        auto superBlockId = idx >> level1_size;
        auto bitId        = idx &  63;
        uint64_t a={};
        for (size_t i{0}; i<= symb; ++i) {
            a +=   level0[level0Id][i]
                 + level1[level1Id][i]
//                 + level2[level2Id][i]
                 + superBlocks[superBlockId][i];

        }
        return bits[level0Id].prefix_rank(bitId, symb) + a;
    }


    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
        prefetch(idx);

        auto level0Id     = idx >>  6;
        auto level1Id     = idx >> level0_size;
//        auto level2Id     = idx >> level1_size;
        auto superBlockId = idx >> level1_size;
        auto bitId        = idx &  63;
        auto res = std::array<uint64_t, TSigma>{};
        for (size_t symb{0}; symb < TSigma; ++symb) {
            res[symb] =   bits[level0Id].rank(bitId, symb)
                        + level0[level0Id][symb]
                        + level1[level1Id][symb]
//                        + level2[level2Id][symb]
                        + superBlocks[superBlockId][symb]
                        + C[symb];
        }
        return res;
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, TSigma>, std::array<uint64_t, TSigma>> {
        prefetch(idx);

        auto level0Id     = idx >>  6;
        auto level1Id     = idx >> level0_size;
//        auto level2Id     = idx >> level1_size;
        auto superBlockId = idx >> level1_size;
        auto bitId        = idx &  63;

        std::array<uint64_t, TSigma> prs;
        auto rs = bits[level0Id].all_ranks(bitId);

        rs[0] +=   level0[level0Id][0]
                 + level1[level1Id][0]
//                 + level2[level2Id][0]
                 + superBlocks[superBlockId][0]
                 + C[0];

        prs[0] = rs[0];
        for (size_t symb{1}; symb < TSigma; ++symb) {
            auto a =   level0[level0Id][symb]
                     + level1[level1Id][symb]
//                     + level2[level2Id][symb]
                     + superBlocks[superBlockId][symb];

            prs[symb] = prs[symb-1] + rs[symb] + a;
            rs[symb] += C[symb] + a;
        }
        return {rs, prs};
    }

    size_t symbol(uint64_t idx) const {
        auto level0Id     = idx >>  6;
        auto bitId        = idx &  63;
        return bits[level0Id].symbol(bitId);
    }

    size_t rank_symbol(uint64_t idx) const {
        prefetch(idx);

        auto level0Id     = idx >>  6;
        auto level1Id     = idx >> level0_size;
//        auto level2Id     = idx >> level1_size;
        auto superBlockId = idx >> level1_size;
        auto bitId = idx & 63;

        auto [rank, symb] = bits[level0Id].rank_symbol(bitId);
        return    rank
                + level0[level0Id][symb]
                + level1[level1Id][symb]
//                + level2[level2Id][symb]
                + superBlocks[superBlockId][symb]
                + C[symb];
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bits, level0, level1, superBlocks, C);
    }
};


template <size_t TSigma, size_t TAlignment>
struct OccTable {
    using TLengthType = uint64_t;
    static constexpr size_t Sigma = TSigma;

    Bitvector<Sigma, TAlignment> bitvector;

    static size_t expectedMemoryUsage(size_t length) {
        using Block = typename Bitvector<TSigma, TAlignment>::BlockL1;
        auto blockSize = std::max(alignof(Block), sizeof(Block));

        size_t C           = sizeof(uint64_t) * (Sigma+1);
        size_t blocks      = blockSize        * (length+1) / 64;
        size_t superblocks = sizeof(uint64_t) * (length+1) / (1ul << (sizeof(uint16_t) * 8));
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

    size_t rank_symbol(size_t idx) const {
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

namespace eprV5 {
template <size_t TSigma>
struct OccTable : eprV5_impl::OccTable<TSigma, 1> {
    static auto name() -> std::string {
        return "EPR V5";
    }

    static auto extension() -> std::string {
        return "eprv5";
    }
};
static_assert(checkOccTable<OccTable>);
}



}
}
