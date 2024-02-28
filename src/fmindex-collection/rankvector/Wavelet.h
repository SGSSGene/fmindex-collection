// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../builtins.h"
#include "../bitvector/Bitvector.h"
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
template <size_t TSigma, BitVector_c Bitvector = ::fmindex_collection::bitvector::Bitvector>
struct Wavelet {
    static constexpr size_t Sigma = TSigma;

    static constexpr auto bits = required_bits(TSigma);
    static constexpr auto bvct = pow(2, bits);

    std::array<Bitvector, bvct> bitvector;
    size_t                      totalLength;

    Wavelet() = default;
    Wavelet(std::span<uint8_t const> _symbols) {
        totalLength = _symbols.size();
        auto length = _symbols.size();

        auto which_bv = [](uint64_t symb, auto cb) {
            uint64_t id{0};
            uint64_t factor{1};
            uint64_t mask = 1u << (bits-1);
            for (uint64_t b{0}; b < bits; ++b) {
                auto bit = (symb & mask) != 0;
                cb(id, bit);
                id += (bit + 1) * factor;
                factor = factor * 2;
                mask = mask >> 1;
            }
        };


        for (uint64_t size{1}; size <= length; ++size) {
            auto symb = _symbols[size-1];
            which_bv(symb, [&](uint64_t id, uint64_t bit) {
                bitvector[id].push_back(bit);
            });
        }
    }

    size_t size() const {
        return totalLength;
    }

    uint8_t symbol(uint64_t idx) const {
        for (uint64_t i{0}; i < Sigma-1; ++i) {
            if (rank(idx+1, i) != rank(idx, i)) {
                return i;
            }
        }
        return Sigma-1;
    }

    uint64_t rank(uint64_t idx, uint8_t symb) const {
        auto which_bv = [](uint64_t symb, auto cb) {
            uint64_t id{0};
            uint64_t factor{1};
            uint64_t mask = 1u << (bits-1);
            for (uint64_t b{0}; b < bits; ++b) {
                auto bit = (symb & mask) != 0;
                cb(id, bit);
                id += (bit + 1) * factor;
                factor = factor * 2;
                mask = mask >> 1;
            }
        };
        uint64_t a = idx;
        which_bv(symb, [&](uint64_t id, uint64_t value) {
            auto newIdx = bitvector[id].rank(a);
            if (value == 0) {
                a = a - newIdx;
            } else {
                a = newIdx;
            }
        });
        return a;
    }

    uint64_t prefix_rank(uint64_t idx, uint8_t symb) const {
        auto which_bv = [](size_t symb, auto cb) {
            uint64_t id{0};
            uint64_t factor{1};
            uint64_t mask = 1u << (bits-1);
            for (uint64_t b{0}; b < bits; ++b) {
                auto bit = (symb & mask) != 0;
                cb(id, bit);
                id += (bit + 1) * factor;
                factor = factor * 2;
                mask = mask >> 1;
            }
        };
        uint64_t a{};
        uint64_t pos = idx;
        which_bv(symb+1, [&](uint64_t id, uint8_t value) {
            auto newIdx = bitvector[id].rank(pos);
            if (value == 0) {
                pos = pos - newIdx;
            } else {
                a += pos - newIdx;
                pos = newIdx;
            }
        });

        return a;
    }

    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, Sigma> {
        std::array<uint64_t, Sigma> rs{0};
        for (uint64_t i{0}; i < Sigma; ++i) {
            rs[i] = rank(idx, i);
        }
        return rs;
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, Sigma>, std::array<uint64_t, Sigma>> {
        std::array<uint64_t, Sigma> rs{0};
        std::array<uint64_t, Sigma> prs{0};
        for (uint64_t i{0}; i < Sigma; ++i) {
            rs[i] = rank(idx, i);
            prs[i] = prefix_rank(idx, i);
        }
        return {rs, prs};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bitvector, totalLength);
    }

};
template <uint64_t TSigma>
using Wavelet_Default = Wavelet<TSigma>;

static_assert(checkRankVector<Wavelet_Default>);

}
