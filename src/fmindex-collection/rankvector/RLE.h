// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../builtins.h"
#include "../bitvector/CompactBitvector.h"
#include "Naive.h"
#include "concepts.h"

#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <span>
#include <vector>

/**
 * A run length encoded symbols, see unpublished paper
 *
 */
namespace fmindex_collection::rankvector {

template <size_t TSigma, size_t encodingBlockSize, typename RankVector = Naive<TSigma>, typename RecRankVector = Naive<TSigma>>
struct RLE {
    using BitVector  = bitvector::CompactBitvector;
    static constexpr size_t Sigma = TSigma;

//    size_t        encodingBlockSize;
    RankVector    bitvector1{};
    RecRankVector bitvector2{};
    BitVector     partition{};

    RLE() = default;
    RLE(std::span<uint8_t const> _symbols/*, size_t _encodingBlockSize*/)
//        : encodingBlockSize{_encodingBlockSize}
    {
        assert(encodingBlockSize > 1);
        auto symbols1     = std::vector<uint8_t>{};
        auto symbols2     = std::vector<uint8_t>{};
        auto partitionSym = std::vector<uint8_t>{};


        size_t trailingChars = (_symbols.size() % encodingBlockSize);
        for (size_t i{0}; i < _symbols.size() - trailingChars; i += encodingBlockSize) {
            bool allTheSame = true;
            for (size_t j{1}; j < encodingBlockSize; ++j) {
                allTheSame &= (_symbols[i] == _symbols[i+j]);
            }
            if (allTheSame) {
                partitionSym.push_back(1);
                symbols2.push_back(_symbols[i]);
            } else {
                partitionSym.push_back(0);
                for (size_t j{0}; j < encodingBlockSize; ++j) {
                    symbols1.push_back(_symbols[i+j]);
                }
            }
        }

        if (trailingChars > 0) {
            partitionSym.push_back(0);
            for (size_t i{_symbols.size()-trailingChars}; i < _symbols.size(); i += 1) {
                symbols1.push_back(_symbols[i]);
            }
        }

        bitvector1 = RankVector{symbols1};
        bitvector2 = RecRankVector{symbols2};
        partition  = BitVector(partitionSym.size(), [&](size_t i) {
            return partitionSym[i];
        });
    }

    size_t size() const {
        return bitvector1.size()
            + bitvector2.size() * encodingBlockSize;
    }

    uint8_t symbol(uint64_t idx) const {
        auto nbr = partition.rank(idx / encodingBlockSize);
        auto s = partition.symbol(idx / encodingBlockSize);
        if (s == 1) {
            return bitvector2.symbol(nbr);
        }
        return bitvector1.symbol(idx - nbr*encodingBlockSize);
    }

    uint64_t rank(uint64_t idx, uint64_t symb) const {
        auto nbr = partition.rank(idx / encodingBlockSize);
        auto v2  = bitvector2.rank(nbr, symb) * encodingBlockSize;
        auto s = partition.symbol(idx/encodingBlockSize);
        if (s == 0) {
            auto v = bitvector1.rank(idx - nbr*encodingBlockSize, symb);
            return v + v2;
        } else {
            auto tail = idx % encodingBlockSize;
            auto v = bitvector1.rank(idx - nbr*encodingBlockSize - tail, symb);
            if (bitvector2.symbol(nbr) == symb) {
                v = v + tail;
            }
            return v + v2;
        }
    }

    uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
        auto nbr = partition.rank(idx / encodingBlockSize);
        auto v2  = bitvector2.prefix_rank(nbr, symb) * encodingBlockSize;

        auto s = partition.symbol(idx/encodingBlockSize);
        if (s == 0) {
            auto v = bitvector1.prefix_rank(idx - nbr*encodingBlockSize, symb);
            return v + v2;
        } else {
            auto tail = idx % encodingBlockSize;
            auto v = bitvector1.prefix_rank(idx - nbr*encodingBlockSize - tail, symb);

            if (bitvector2.symbol(nbr) <= symb) {
                v = v + tail;
            }
            return v + v2;
        }
    }

    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
        auto nbr = partition.rank(idx / encodingBlockSize);
        auto v2  = bitvector2.all_ranks(nbr);
        auto s = partition.symbol(idx/encodingBlockSize);
        if (s == 0) {
            auto v = bitvector1.all_ranks(idx - nbr*encodingBlockSize);
            for (size_t i{0}; i < TSigma; ++i) {
                v2[i] = v[i] + v2[i] * encodingBlockSize;
            }
        } else {
            auto tail = idx % encodingBlockSize;
            auto v = bitvector1.all_ranks(idx - nbr*encodingBlockSize - tail);
            for (size_t i{0}; i < TSigma; ++i) {
                v2[i] = v[i] + v2[i] * encodingBlockSize;
            }
            auto curSymb = bitvector2.symbol(nbr);
            v2[curSymb] += tail;
        }
        return v2;
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, TSigma>, std::array<uint64_t, TSigma>> {
        auto rs = all_ranks(idx);
        auto prs = rs;
        for (size_t i{1}; i < TSigma; ++i) {
            prs[i] = prs[i-1] + prs[i];
        }
        return {rs, prs};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(encodingBlockSize, bitvector1, bitvector2, partition);
    }

};

template <size_t TSigma> using RLEInstance = RLE<TSigma, 4>;
static_assert(checkRankVector<RLEInstance>);

template <uint64_t TSigma, size_t encodingBlockSize, typename RankVector = Naive<TSigma>, size_t depth = 0>
struct rRLE : RLE<TSigma, encodingBlockSize, RankVector, rRLE<TSigma, encodingBlockSize, RankVector, depth-1>>
{};

template <uint64_t TSigma, size_t encodingBlockSize, typename RankVector>
struct rRLE<TSigma, encodingBlockSize, RankVector, 0> : RLE<TSigma, encodingBlockSize, RankVector, RankVector>
{};

template <size_t TSigma> using rRLEInstance = rRLE<TSigma, 2, Naive<TSigma>, 2>;
static_assert(checkRankVector<rRLEInstance>);

}
