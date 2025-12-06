// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../builtins.h"
#include "../bitvector/Bitvector.h"
#include "../bitvector/PrunedBitvector.h"
#include "concepts.h"

#include <algorithm>
#include <array>
#include <bit>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

namespace fmc::string {

/* Implements the concept `String_c`
 *
 * \param TSigma size of the alphabet
 */
template <size_t TSigma, Bitvector_c Bitvector = ::fmc::bitvector::Bitvector>
struct Wavelet {
    static constexpr size_t Sigma = TSigma;

    static constexpr auto bits = std::bit_width(TSigma-1);
    static constexpr auto bvct = std::bit_ceil(TSigma);

    std::array<Bitvector, bvct> bitvector;
    size_t                      totalLength;

    static auto extractBitsFromSymb(uint8_t symb) -> std::array<std::tuple<uint8_t, uint8_t>, bits> {
        auto res = std::array<std::tuple<uint8_t, uint8_t>, bits>{};

        for (uint8_t b{0}; b < bits; ++b) {
            auto bitId          = bits - b - 1;
            uint8_t bit         = (symb >> bitId) & 1;
            uint8_t id_offset   = (1ull<<b)-1;
            uint8_t symb_offset = symb >> (bitId+1);
            uint8_t id = id_offset + symb_offset;
            res[b] = {bit, id};

        }
        return res;
    }
    std::array<std::array<std::tuple<uint8_t, uint8_t>, bits>, Sigma> lut{[]() {
        auto lut = std::array<std::array<std::tuple<uint8_t, uint8_t>, bits>, Sigma>{};
        for (size_t s{0}; s < Sigma; ++s) {
            lut[s] = extractBitsFromSymb(s);
        }
        return lut;
    }()};


    Wavelet() = default;
    Wavelet(std::span<uint8_t const> _symbols)
        : totalLength {_symbols.size()} {
        for (uint64_t size{0}; size < totalLength; ++size) {
            auto symb = _symbols[size];
            auto const& res = lut[symb];
//            auto res = extractBitsFromSymb(symb);
            for (auto [bit, id] : res) {
                bitvector[id].push_back(bit);
            }
        }
    }

    size_t size() const {
        return totalLength;
    }

    uint8_t symbol(uint64_t idx) const {
        assert(idx < totalLength);
        uint8_t symb{};
        for (uint8_t b{0}; b < bits; ++b) {
            uint8_t id_offset = (1ull<<b)-1;
            uint8_t id        = id_offset + symb;

            auto [bit, newIdx] = [&]() -> std::tuple<size_t, size_t> {
                if (id < bitvector.size()) {
                    size_t bit    = bitvector[id].symbol(idx);
                    size_t newIdx = bitvector[id].rank(idx);
                    return {bit, newIdx};
                } else {
                    return {0, 0};
                }
            }();
            symb = (symb<<1) | bit;
            if (bit == 0) {
                idx = idx - newIdx;
            } else {
                idx = newIdx;
            }

        }
        return symb;
    }

    uint64_t rank(uint64_t idx, uint8_t symb) const {
        assert(idx <= totalLength);
        assert(symb < TSigma);

        auto const& res = lut[symb];
//        auto res = extractBitsFromSymb(symb);
        for (auto [bit, id] : res) {
            size_t newIdx = bitvector[id].rank(idx);
            if (bit == 0) {
                idx = idx - newIdx;
            } else {
                idx = newIdx;
            }
        }
        return idx;
    }

    uint64_t prefix_rank(uint64_t idx, uint8_t symb) const {
        assert(idx <= totalLength);
        assert(symb <= TSigma);

        if (symb == 0) return 0;
        symb -= 1;

        uint64_t a{};
        auto const& res = lut[symb];
//        auto res = extractBitsFromSymb(symb);
        for (auto [bit, id] : res) {
            size_t newIdx = bitvector[id].rank(idx);
            if (bit == 0) {
                idx = idx - newIdx;
            } else {
                a = a + idx - newIdx;
                idx = newIdx;
            }
        }
        return a + idx;
    }

    template <size_t b=0, size_t symb=0, typename CB>
    void bv_all(size_t idx, CB const& cb) const {
        if constexpr(symb >= Sigma) {
            return;
        } else if constexpr (b == bits) {
            cb.template operator()<symb>(idx);
        } else {
            uint8_t id_offset = (1ull<<b)-1;
            uint8_t id        = id_offset + symb;

            auto newIdx = [&]() -> size_t {
                return bitvector[id].rank(idx);
            }();


            bv_all<b+1, symb<<1>(idx-newIdx, cb);
            bv_all<b+1, (symb<<1) | 1>(newIdx, cb);
        }
    }


    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, Sigma> {
        assert(idx <= totalLength);

        auto rs = std::array<uint64_t, Sigma>{};

        bv_all(idx, [&]<size_t symb>(size_t count) {
            rs[symb] = count;
        });
        return rs;
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, Sigma>, std::array<uint64_t, Sigma>> {
        assert(idx <= totalLength);

        auto rs  = all_ranks(idx);
        auto prs = std::array<uint64_t, Sigma>{};
        for (uint64_t i{1}; i < Sigma; ++i) {
            prs[i] = prs[i-1] + rs[i-1];
        }
        return {rs, prs};
    }

    template <typename Archive>
    void serialize(this auto&& self, Archive& ar) {
        ar(self.bitvector, self.totalLength);
    }
};

template <uint64_t TSigma>
using Wavelet_Default = Wavelet<TSigma>;
static_assert(checkString_c<Wavelet_Default>);

template <size_t Sigma>
struct PrunedWavelet : Wavelet<Sigma, bitvector::PrunedBitvector> {
    using Wavelet<Sigma, bitvector::PrunedBitvector>::Wavelet;
};

template <uint64_t TSigma>
using PrunedWavelet_Default = PrunedWavelet<TSigma>;

static_assert(checkString_c<PrunedWavelet_Default>);


}
