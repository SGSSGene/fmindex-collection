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

template <size_t TSigma, size_t encodingBlockSize, size_t DifferentBlocks = 16, template <size_t, typename...> typename String = Wavelet, template <size_t, typename...> typename MixedString = String>
struct RBBwtV2 {
    static constexpr size_t Sigma = TSigma;

    String<DifferentBlocks> topLevelVector{};
    MixedString<TSigma>     mixedLevelVector{};

    std::array<std::array<size_t, DifferentBlocks>, TSigma> factors{};
    std::array<std::array<uint8_t, encodingBlockSize>, DifferentBlocks> prefix{};
    std::array<std::array<std::array<uint8_t, TSigma>, encodingBlockSize+1>, DifferentBlocks> prefix2{};

    size_t totalSize{};

    RBBwtV2() = default;
    RBBwtV2(std::span<uint8_t const> _symbols) {
        auto blockName = [&](size_t offset) {
            size_t a{};
            for (size_t i{0}; i < encodingBlockSize; ++i) {
                a = a * TSigma + _symbols[offset+encodingBlockSize-i-1];
            }
            return a;
        };

        // collects all blocks and make a list of topBlocks
        auto topBlocks = std::unordered_map<size_t, size_t>{};
        {
            auto map = std::unordered_map<size_t, size_t>{};
            for (size_t i{0}; (i+encodingBlockSize-1) < _symbols.size(); i += encodingBlockSize) {
                map[blockName(i)] += 1;
            }
            auto sortedBlocks = std::vector<std::tuple<size_t, size_t>>{};
            for (auto const& [id, count] : map) {
                sortedBlocks.emplace_back(count, id);
            }
            std::ranges::sort(sortedBlocks);
            std::ranges::reverse(sortedBlocks);
            for (size_t i{0}; i < sortedBlocks.size() && i < DifferentBlocks-1; ++i) {
                auto [count, id] = sortedBlocks[i];
                topBlocks[id] = i;
                // create factors
                auto acc = std::array<uint8_t, TSigma>{};
                for (size_t j{0}; j < encodingBlockSize; ++j) {
                    auto v = id % TSigma;
                    factors[v][i] += 1;
                    prefix2[i][j] = acc;
                    acc[v] += 1;
                    prefix[i][j] = v;
                    id = id / TSigma;
                }
                prefix2[i][encodingBlockSize] = acc;
            }
        }
        // naming of blocks
        {
            auto topLevelNames = std::vector<uint8_t>{};
            auto mixedLevel    = std::vector<uint8_t>{};

            size_t trailingChars = (_symbols.size() % encodingBlockSize);
            for (size_t i{0}; i < _symbols.size() - trailingChars; i += encodingBlockSize) {
                auto name = blockName(i);
                if (auto iter = topBlocks.find(name); iter != topBlocks.end()) {
                    topLevelNames.emplace_back(iter->second);
                } else {
                    topLevelNames.emplace_back(DifferentBlocks-1);
                    for (size_t j{0}; j < encodingBlockSize; ++j) {
                        mixedLevel.push_back(_symbols[i+j]);
                    }
                }
            }

            if (trailingChars > 0) {
                topLevelNames.emplace_back(DifferentBlocks-1);
                for (size_t i{_symbols.size() - trailingChars}; i < _symbols.size(); ++i) {
                    mixedLevel.push_back(_symbols[i]);
                }
            }

            topLevelVector   = String<DifferentBlocks>{topLevelNames};
            mixedLevelVector = MixedString<TSigma>{mixedLevel};

/*            for (size_t i{0}; i < topLevelNames.size(); ++i) {
                assert(topLevelVector.symbol(i) == topLevelNames[i]);
            }
            {
                auto acc = std::array<size_t, DifferentBlocks>{};
                for (size_t i{0}; i  < topLevelNames.size(); ++i) {
                    std::cout << i << " ";
                    for (size_t s{}; s < DifferentBlocks; ++s) {
                        auto v = topLevelVector.rank(i, s);
                        std::cout << v << " ";
                        assert(v == acc[s]);
                        acc[s] += (topLevelNames[i] == s);
                    }
                    std::cout << "\n";
                }
            }*/
        }
        totalSize = _symbols.size();
    }

    size_t size() const {
        return totalSize;
    }

    uint8_t symbol(uint64_t idx) const {
        assert(idx < totalSize);

        auto indicator = topLevelVector.symbol(idx/encodingBlockSize);
        if (indicator == DifferentBlocks-1) { // mixed Level
            auto nbr = topLevelVector.rank(idx/encodingBlockSize, DifferentBlocks-1);
            return mixedLevelVector.symbol((nbr*encodingBlockSize) + (idx%encodingBlockSize));
        } else {               // inside a 'top' block
            return prefix[indicator][idx%encodingBlockSize];
        }
    }

    uint64_t rank(uint64_t idx, uint64_t symb) const {
        assert(idx <= totalSize);
        assert(symb < Sigma);

        if (idx == 0) return 0;
//        for (size_t i{0}; i  < totalSize / encodingBlockSize; ++i) {
//            std::cout << i << " ";
//            for (size_t s{}; s < DifferentBlocks; ++s) {
//                auto v = topLevelVector.rank(i, s);
//                std::cout << v << " ";
//            }
//            std::cout << "\n";
//        }

        size_t r{};
        auto indicators = topLevelVector.all_ranks(idx/encodingBlockSize);
        auto blockType  = topLevelVector.symbol((idx-1)/encodingBlockSize);

        auto nbr = indicators[DifferentBlocks-1];

        size_t prefixLen = idx % encodingBlockSize;
        r +=  mixedLevelVector.rank((nbr*encodingBlockSize) + ((blockType == DifferentBlocks-1)?prefixLen:0), symb);

        if (blockType != DifferentBlocks-1) {
            r += prefix2[blockType][prefixLen][symb];
        }
        for (size_t i{0}; i < DifferentBlocks-1; ++i) {
            r += factors[symb][i] * indicators[i];
        }

        return r;
    }

    uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
        assert(idx <= totalSize);
        assert(symb <= Sigma);

        if (idx == 0) return 0;
        size_t r{};
        auto indicators = topLevelVector.all_ranks(idx/encodingBlockSize);
        auto blockType  = topLevelVector.symbol((idx-1)/encodingBlockSize);

        auto nbr = indicators[DifferentBlocks-1];

        size_t prefixLen = idx % encodingBlockSize;
        r +=  mixedLevelVector.prefix_rank((nbr*encodingBlockSize) + ((blockType == DifferentBlocks-1)?prefixLen:0), symb);

        if (blockType != DifferentBlocks-1) {
            for (size_t s{0}; s <= symb; ++s) {
                r += prefix2[blockType][prefixLen][s];
            }
        }
        for (size_t s{0}; s <= symb; ++s) {
            for (size_t i{0}; i < DifferentBlocks-1; ++i) {
                r += factors[s][i] * indicators[i];
            }
        }
        return r;
    }

    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
        assert(idx <= totalSize);

        auto res = std::array<uint64_t, TSigma>{};
        if (idx == 0) return res;
        for (size_t s{0}; s < TSigma; ++s) {
            res[s] = rank(idx, s);
        }
        return res;
/*        auto nbr = partition.rank(idx / encodingBlockSize);
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
        return v2;*/
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, TSigma>, std::array<uint64_t, TSigma>> {
        assert(idx <= totalSize);

        auto rs = all_ranks(idx);
        auto prs = std::array<uint64_t, TSigma>{};
        for (size_t i{1}; i < TSigma; ++i) {
            prs[i] = prs[i-1] + rs[i-1];
        }
        return {rs, prs};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
//        ar(topLevelVector);
//        ar(mixedLevelVector);
        ar(topLevelVector, mixedLevelVector, factors, prefix, prefix2, totalSize);
    }

};

template <size_t TSigma> using RBBwtV2Instance = RBBwtV2<TSigma, 4>;
static_assert(checkString_c<RBBwtV2Instance>);

}
