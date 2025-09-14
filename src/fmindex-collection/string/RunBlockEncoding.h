// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../builtins.h"
#include "../bitvector/CompactBitvector.h"
#include "Wavelet.h"
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
namespace fmc::string {

template <size_t TSigma, size_t encodingBlockSize, String_c String = Wavelet<TSigma>, String_c RecString = Wavelet<TSigma>>
struct RunBlockEncoding {
    static_assert(String::Sigma == RecString::Sigma);

    using BitVector = bitvector::CompactBitvector;
    static constexpr size_t Sigma = TSigma;

//    size_t        encodingBlockSize;
    String     bitvector1{};
    RecString  bitvector2{};
    BitVector  partition{};

    RunBlockEncoding() = default;
    RunBlockEncoding(std::span<uint8_t const> _symbols/*, size_t _encodingBlockSize*/)
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
        symbols1.push_back(0);
        symbols2.push_back(0);

        bitvector1 = String{symbols1};
        bitvector2 = RecString{symbols2};
        partition  = BitString(partitionSym.size(), [&](size_t i) {
            return partitionSym[i];
        });
    }

    size_t size() const {
        return bitvector1.size() - 1
            + (bitvector2.size() - 1) * encodingBlockSize;
    }

    uint8_t symbol(uint64_t idx) const {
        assert(idx < size());

        auto nbr = partition.rank(idx / encodingBlockSize);
        auto s = partition.symbol(idx / encodingBlockSize);
        if (s == 1) {
            return bitvector2.symbol(nbr);
        }
        return bitvector1.symbol(idx - nbr*encodingBlockSize);
    }

    uint64_t rank(uint64_t idx, uint64_t symb) const {
        assert(idx <= size());
        assert(symb < Sigma);

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
        assert(idx <= size());
        assert(symb <= Sigma);

        auto nbr = partition.rank(idx / encodingBlockSize);
        auto v2  = bitvector2.prefix_rank(nbr, symb) * encodingBlockSize;

        auto s = partition.symbol(idx/encodingBlockSize);
        if (s == 0) {
            auto v = bitvector1.prefix_rank(idx - nbr*encodingBlockSize, symb);
            return v + v2;
        } else {
            auto tail = idx % encodingBlockSize;
            auto v = bitvector1.prefix_rank(idx - nbr*encodingBlockSize - tail, symb);

            if (bitvector2.symbol(nbr) < symb) {
                v = v + tail;
            }
            return v + v2;
        }
    }

    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
        assert(idx <= size());

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
        assert(idx <= size());

        auto rs = all_ranks(idx);
        auto prs = std::array<uint64_t, TSigma>{};
        for (size_t i{1}; i < TSigma; ++i) {
            prs[i] = prs[i-1] + rs[i-1];
        }
        return {rs, prs};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(encodingBlockSize, bitvector1, bitvector2, partition);
    }

};

template <size_t TSigma> using RunBlockEncodingInstance = RunBlockEncoding<TSigma, 4>;
static_assert(checkString_c<RunBlockEncodingInstance>);

template <uint64_t TSigma, size_t encodingBlockSize, String_c String = Wavelet<TSigma>, size_t depth = 0>
struct rRunBlockEncoding : RunBlockEncoding<TSigma, encodingBlockSize, String, rRunBlockEncoding<TSigma, encodingBlockSize, String, depth-1>>
{};

template <uint64_t TSigma, size_t encodingBlockSize, String_c String>
struct rRunBlockEncoding<TSigma, encodingBlockSize, String, 0> : RunBlockEncoding<TSigma, encodingBlockSize, String, String>
{};

template <size_t TSigma> using rRunBlockEncodingInstance = rRunBlockEncoding<TSigma, 2, Wavelet<TSigma>, 2>;
static_assert(checkString_c<rRunBlockEncodingInstance>);

}
