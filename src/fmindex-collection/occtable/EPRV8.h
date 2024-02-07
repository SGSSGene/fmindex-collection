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
#include <cmath>
#include <cstdint>
#include <span>
#include <vector>


// Same as EPRV5 but with runtime sigma instead of compile time


namespace fmindex_collection::occtable {
namespace eprV8_impl {

struct Bitvector {

    size_t sigma{};
    size_t bitct{size_t(std::ceil(std::log2(sigma)))};
    size_t bvct{size_t(std::exp2(bitct))};

    // number of full length bitvectors needed `2^bitct ≥ TSigma`
//    static constexpr auto bitct = required_bits(TSigma-1);
    // next full power of 2
//    static constexpr auto bvct  = pow(2, bitct);

    struct InBitsView {
        uint64_t const* bits;
        size_t const bitct;
//        std::array<uint64_t, bitct> bits{};

        uint64_t rank(uint64_t idx, uint64_t symb) const {
            assert(idx < 64);

            size_t mask = -1;
            for (size_t i{0}; i < bitct; ++i) {
                mask &= bits[i] ^ -((~symb>>i)&1);
            }

            auto bitset = std::bitset<64>(mask) << (64-idx);
            return bitset.count();
        }

        uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
            assert(idx < 64);

            size_t mask = 0;
            for (size_t _symb{0}; _symb <= symb; ++_symb) {
                size_t imask = -1;
                for (size_t i{0}; i < bitct; ++i) {
                    imask &= bits[i] ^ -((~symb>>i)&1);
                }
                mask |= imask;
            };

            auto bitset = std::bitset<64>(mask) << (64-idx);
            return bitset.count();
        }

        void all_ranks(uint64_t idx, std::span<uint64_t> out_rs) const {
            assert(idx < 64);

            std::ranges::fill(out_rs, 0);

            for (uint64_t i{0}; i < out_rs.size(); ++i) {
                size_t mask = -1;
                for (uint64_t j{0}; j < bitct; ++j) {
                    mask &= bits[j] ^ -((~i>>j)&1);
                }
                out_rs[i] = (std::bitset<64>(mask) << (64 - idx)).count();
            }
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

            for (size_t i{0}; i < bitct; ++i) {
                auto b = (bits[i] >> idx) & 1;
                mask |= bits[i] ^ -b;
                symb |= b << i;
            }

            auto bitset = std::bitset<64>{~mask} << (64-idx);

            return {bitset.count(), symb};
        }
    };

    using blockL0_t = uint8_t;
    using blockL1_t = uint16_t;
    static constexpr uint64_t level0_size = sizeof(blockL0_t) * 8;
    static constexpr uint64_t level1_size = sizeof(blockL1_t) * 8;

    using BlockL0 = blockL0_t;
    using BlockL1 = blockL1_t;

    std::vector<uint64_t> bits;
//    std::vector<InBitsView> bitsView;
//    std::vector<InBits> bits;
    std::vector<BlockL0> level0;
    std::vector<BlockL1> level1;
    std::vector<uint64_t> superBlocks;


    std::vector<uint64_t> C;
//    std::array<uint64_t, TSigma+1> C;

    Bitvector() = default;
    Bitvector(std::span<uint8_t const> _bwt, size_t _sigma)
        : sigma{_sigma}
    {
        auto const length = _bwt.size();
        level1.reserve(length*sigma/(1ull<<level1_size)+2);
        level0.reserve(length*sigma/64+2);
        bits.reserve(length*bitct/64+2);

        std::vector<blockL0_t> blockL0_acc(sigma, 0);
        std::vector<blockL1_t> blockL1_acc(sigma, 0);
        std::vector<uint64_t> sblock_acc(sigma, 0);


        for (uint64_t size{0}; size < length+1; ++size) {
            if (size % (1ull<<level1_size) == 0) { // new l3 block
                for (auto a : sblock_acc) {
                    superBlocks.emplace_back(a);
                }
                for (size_t i{0}; i < sigma; ++i) {
                    level1.emplace_back();
                    level0.emplace_back();
                }
                for (size_t i{0}; i < bitct; ++i) bits.emplace_back();
//                bits.emplace_back();

                std::ranges::fill(blockL0_acc, 0);
                std::ranges::fill(blockL1_acc, 0);
//                blockL0_acc = {};
//                blockL1_acc = {};
            } else if (size % (1ull<<level0_size) == 0) { // new l1 block
                for (auto a : blockL1_acc) {
                    level1.emplace_back(a);
                }
                for (size_t i{0}; i < sigma; ++i) level0.emplace_back();
                for (size_t i{0}; i < bitct; ++i) bits.emplace_back();
//                bits.emplace_back();

                std::ranges::fill(blockL0_acc, 0);
//                blockL0_acc = {};
            } else if (size % 64 == 0) { // new l0 block
                for (auto a : blockL0_acc) {
                    level0.emplace_back(a);
                }
                for (size_t i{0}; i < bitct; ++i) bits.emplace_back();
//                bits.emplace_back();
//                level0.back() = blockL0_acc;
            }
            // Abort, we only wanted to add a new block to the end, if required
            if (size == length) continue;

            auto level0Id     = size >>  6;
            auto bitId        = size &  63;

            uint64_t symb = _bwt[size];

            for (uint64_t i{}; i < bitct; ++i) {
                auto b = ((symb>>i)&1);
                bits[level0Id*bitct+i] |= (b << bitId);
            }
            blockL0_acc[symb] += 1;
            blockL1_acc[symb] += 1;
            sblock_acc[symb] += 1;
        }

        C.resize(sigma+1);
        C[0] = 0;
        for (uint64_t i{0}; i < sigma; ++i) {
            C[i+1] = sblock_acc[i] + C[i];
        }
    }

    uint64_t memoryUsage() const {
        return    bits.size() * sizeof(bits.back())
                + level0.size() * sizeof(level0.back())
                + level1.size() * sizeof(level1.back())
                + superBlocks.size() * sizeof(superBlocks.back())
                + sizeof(C);
    }

    void prefetch(uint64_t idx) const {
        auto level0Id     = idx >>  6;
        auto level1Id     = idx >> level0_size;
        auto superBlockId = idx >> level1_size;

        __builtin_prefetch(reinterpret_cast<void const*>(&level1[level1Id]), 0, 0);
        __builtin_prefetch(reinterpret_cast<void const*>(&level0[level0Id]), 0, 0);
        __builtin_prefetch(reinterpret_cast<void const*>(&bits[level1Id]), 0, 0);
        __builtin_prefetch(reinterpret_cast<void const*>(&superBlocks[superBlockId]), 0, 0);
    }

    uint64_t rank(uint64_t idx, uint64_t symb) const {
        prefetch(idx);

        auto level0Id     = idx >>  6;
        auto level1Id     = idx >> level0_size;
        auto superBlockId = idx >> level1_size;
        auto bitId        = idx &  63;
        return    InBitsView{&bits[level0Id*bitct], bitct}.rank(bitId, symb)
                + level0[level0Id*sigma+symb]
                + level1[level1Id*sigma+symb]
                + superBlocks[superBlockId*sigma+symb]
                + C[symb];
    }

    uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
        prefetch(idx);

        auto level0Id     = idx >>  6;
        auto level1Id     = idx >> level0_size;
        auto superBlockId = idx >> level1_size;
        auto bitId        = idx &  63;
        uint64_t a={};
        for (uint64_t i{0}; i<= symb; ++i) {
            a +=   level0[level0Id*sigma+i]
                 + level1[level1Id*sigma+i]
                 + superBlocks[superBlockId*sigma+i];

        }
        return InBitsView{&bits[level0Id*bitct], bitct}.prefix_rank(bitId, symb) + a;
    }


    auto all_ranks(uint64_t idx, std::span<uint64_t> out) const {
        prefetch(idx);

        auto level0Id     = idx >>  6;
        auto level1Id     = idx >> level0_size;
        auto superBlockId = idx >> level1_size;
        auto bitId        = idx &  63;
        for (uint64_t symb{0}; symb < out.size(); ++symb) {
            out[symb] =   InBitsView{&bits[level0Id*bitct], bitct}.rank(bitId, symb)
                        + level0[level0Id*sigma+symb]
                        + level1[level1Id*sigma+symb]
                        + superBlocks[superBlockId*sigma+symb]
                        + C[symb];
        }
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx, std::span<uint64_t> out_rs, std::span<uint64_t> out_prs) const {
        prefetch(idx);

        auto level0Id     = idx >>  6;
        auto level1Id     = idx >> level0_size;
        auto superBlockId = idx >> level1_size;
        auto bitId        = idx &  63;

        InBitsView{&bits[level0Id*bitct], bitct}.all_ranks(bitId, out_rs);

        out_rs[0] +=   level0[level0Id*sigma]
                     + level1[level1Id*sigma]
                     + superBlocks[superBlockId*sigma]
                     + C[0];

        out_prs[0] = out_rs[0];
        for (uint64_t symb{1}; symb < out_rs.size(); ++symb) {
            auto a =   level0[level0Id*sigma+symb]
                     + level1[level1Id*sigma+symb]
                     + superBlocks[superBlockId*sigma+symb];

            out_prs[symb] = out_prs[symb-1] + out_rs[symb] + a;
            out_rs[symb] += C[symb] + a;
        }
    }

    uint64_t symbol(uint64_t idx) const {
        auto level0Id     = idx >>  6;
        auto bitId        = idx &  63;
        return InBitsView{&bits[level0Id*bitct], bitct}.symbol(bitId);
    }

    uint64_t rank_symbol(uint64_t idx) const {
        prefetch(idx);

        auto level0Id     = idx >>  6;
        auto level1Id     = idx >> level0_size;
        auto superBlockId = idx >> level1_size;
        auto bitId = idx & 63;

        auto [rank, symb] = InBitsView{&bits[level0Id*bitct], bitct}.rank_symbol(bitId);
        return    rank
                + level0[level0Id*sigma+symb]
                + level1[level1Id*sigma+symb]
                + superBlocks[superBlockId*sigma+symb]
                + C[symb];
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bits, level0, level1, superBlocks, C);
        ar(sigma, bitct, bvct);
    }
};


template <size_t TSigma>
struct OccTable {
    static constexpr size_t Sigma = TSigma;

    Bitvector bitvector;

    static uint64_t expectedMemoryUsage(uint64_t length) {
        (void)length;
        return 0;
    }

    OccTable() = default;
    OccTable(std::span<uint8_t const> _bwt)
        : bitvector{_bwt, Sigma}
    {}

    uint64_t memoryUsage() const {
        return bitvector.memoryUsage() + sizeof(OccTable);
    }

    size_t size() const {
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
        auto res = std::tuple<std::array<uint64_t, Sigma>, std::array<uint64_t, Sigma>>{};
        bitvector.all_ranks_and_prefix_ranks(idx, std::get<0>(res), std::get<1>(res));
        return res;
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bitvector);
    }
};


}

namespace eprV8 {
template <size_t TSigma>
struct OccTable : eprV8_impl::OccTable<TSigma> {
    using eprV8_impl::OccTable<TSigma>::OccTable;
    static auto name() -> std::string {
        return "EPR V8";
    }

    static auto extension() -> std::string {
        return "eprv8";
    }
};
static_assert(checkOccTable<OccTable>);
}



}
