// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../builtins.h"
#include "../bitvector/Bitvector.h"
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

namespace fmindex_collection::rankvector {

/* Implements the concept `RankVector`
 *
 * \param TSigma size of the alphabet
 */
template <size_t TSigma, BitVector_c Bitvector = ::fmindex_collection::bitvector::Bitvector>
struct Wavelet2 {
    static constexpr size_t Sigma = TSigma;

    static constexpr auto bits = std::bit_width(TSigma-1);
    static constexpr auto bvct = std::bit_ceil(TSigma);

    std::array<Bitvector, bvct> bitvector;
    size_t                      totalLength;


    struct Node {
        std::unique_ptr<Node> left, right;
//        size_t count{};
        size_t bitvectorId{};

        size_t start;
        size_t len;
        size_t threshold;

        Node(size_t _start, size_t _len, size_t& lastBitvectorId) {
            start = _start;
            len   = _len;
            assert(len > 1);
            bitvectorId = ++lastBitvectorId ;
            if (len == 2) {
                threshold = _start + 1;
//                count = 1;
            } else if (len == 3) {
                threshold = _start + 2;
                left = std::make_unique<Node>(start, 2, lastBitvectorId);
//                count = left->count + 1;
            } else {
                threshold = start + len/2 + len%2;
                auto leftStart  = start;
                auto leftLen    = threshold-start;
                assert(leftLen > 1);
                auto rightStart = threshold;
                auto rightLen   = len - leftLen;
                assert(rightLen > 1);
                left  = std::make_unique<Node>(leftStart, leftLen, lastBitvectorId);
                right = std::make_unique<Node>(rightStart, rightLen, lastBitvectorId);
//                count = left->count + right->count + 1;
            }
        }

        template <typename CB>
        void visit(size_t _symb, CB const& cb) {
            bool bit = _symb >= threshold;
            cb(bitvectorId, bit);
            if (!bit) {
                if (left) {
                    left->visit(_symb, cb);
                }
            } else {
                if (right) {
                    right->visit(_symb, cb);
                }
            }
        }

        template <typename CB1, typename CB2>
        void visitAll(size_t idx, CB1 const& cb1, CB2 const& cb2) const {
            auto rank = cb1(bitvectorId, idx);
            if (!left)  cb2(start, idx-rank);
            else        left->visitAll(idx-rank, cb1, cb2);
            if (!right) cb2(threshold, rank);
            else        right->visitAll(rank, cb1, cb2);
        }


        template <typename CB>
        size_t findSymbol(CB const& cb) const {
            auto bit = cb(bitvectorId);
            if (bit == 0) {
                if (!left) {
                    assert(!right);
                    return start;
                }
                return left->findSymbol(cb);
            } else {
                if (!right) return start+1;
                return right->findSymbol(cb);
            }
        }
    };

    size_t bitvectorCount{};
    Node node{0, Sigma, bitvectorCount};



    auto extractBitsFromSymb(uint8_t symb) const -> std::array<std::tuple<uint8_t, uint8_t>, bits> {
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

/*    void divide(size_t nbr) const -> size_t {
        if (nbr == 1) return 0;
        if (nbr == 2) return 1;

        auto bits = std::bit_width(nbr-1);
        auto leftNbr = pow(2, bits-1);
        auto rightNbr = nbr-leftNbr;
        auto cl = divide(leftNbr);
        auto c2 = divide(rightNbr);
        return c1 + c2 + 1;
    }*/

    template <size_t Sigma, size_t Depth=1>
    auto divide2(size_t symb) const -> std::vector<std::tuple<uint8_t, uint8_t>> {
        if constexpr (Sigma == 1) {
            return {};
        } else if constexpr (Sigma == 2) {
            auto res = std::vector<std::tuple<uint8_t, uint8_t>>{};
            auto bit = (symb >> (bits-Depth)) & 1;
            res.emplace_back(bit, 255); //!TODO
            return res;
        } else {
            auto res = std::vector<std::tuple<uint8_t, uint8_t>>{};
            auto bit = (symb >> (bits-Depth)) & 1;

            constexpr auto bits = std::bit_width(Sigma-1);
            constexpr auto leftNbr = std::bit_ceil(Sigma-1);
            constexpr auto rightNbr = Sigma - leftNbr;
            if (symb < leftNbr) {
                assert(bit == 0);
                res.emplace_back(bit, 255); //!TODO
                for (auto e : divide2<leftNbr, Depth+1>(symb)) {
                    res.push_back(e);
                }
            } else {
                assert(bit == 1);
                res.emplace_back(bit, 255); //!TODO
                for (auto e : divide2<rightNbr, Depth+1>(symb-leftNbr)) {
                    res.push_back(e);
                }
            }

            return res;
        }
    }

    std::array<std::vector<std::tuple<uint8_t, uint8_t>>, Sigma> lut;


    Wavelet2() = default;
    Wavelet2(std::span<uint8_t const> _symbols)
        : totalLength {_symbols.size()} {
        // create lookup
        for (size_t s{0}; s < Sigma; ++s) {
            node.visit(s, [&](size_t id, size_t bit) {
                lut[s].emplace_back(bit, id);
            });
        }

/*        for (size_t i{0}; i < Sigma; ++i) {
            lut[i] = divide2<Sigma>(i);
            auto res = extractBitsFromSymb(i);
            for (size_t j{0}; j < lut[i].size(); ++j) {
                std::get<1>(lut[i][j]) = std::get<1>(res[j]);
                assert(std::get<1>(res[j]) < bitvector.size());
            }
        }*/

        for (uint64_t size{0}; size < totalLength; ++size) {
            auto symb = _symbols[size];
            assert(symb < lut.size());
            auto const& res = lut[symb];
            for (auto [bit, id] : res) {
                assert(id < bitvector.size());
                bitvector[id].push_back(bit);
            }
        }
    }

    size_t size() const {
        return totalLength;
    }

    uint8_t symbol(uint64_t idx) const {
        auto s = node.findSymbol([&](size_t id) {
            auto bit = bitvector[id].symbol(idx);
            auto r   = bitvector[id].rank(idx);

            if (bit == 0) {
                idx = idx - r;
            } else {
                idx = r;
            }

            return bit;
        });
        return s;
/*        uint8_t symb{};
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
                assert(idx >= newIdx);
                idx = idx - newIdx;
            } else {
                idx = newIdx;
            }

        }
        return symb;*/
    }

    uint64_t rank(uint64_t idx, uint8_t symb) const {
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
                if (id < bitvector.size()) {
                    return bitvector[id].rank(idx);
                }
                return 0;
            }();


            bv_all<b+1, symb<<1>(idx-newIdx, cb);
            bv_all<b+1, (symb<<1) | 1>(newIdx, cb);
        }
    }


    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, Sigma> {
        std::array<uint64_t, Sigma> rs{0};
        node.visitAll(idx, [&](size_t id, size_t idx) {
            return bitvector[id].rank(idx);
        }, [&](size_t symb, size_t idx) {
            rs[symb] = idx;
        });
/*        for (size_t s{0}; s < Sigma; ++s) {
            rs[s] = rank(idx, s);
        }*/

/*        bv_all(idx, [&]<size_t symb>(size_t count) {
            rs[symb] = count;
        });*/
        return rs;
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, Sigma>, std::array<uint64_t, Sigma>> {
        auto rs  = all_ranks(idx);
        auto prs = rs;
        for (uint64_t i{1}; i < Sigma; ++i) {
            prs[i] = prs[i-1] + prs[i];
        }
        return {rs, prs};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bitvector, totalLength);
    }

};
template <uint64_t TSigma>
using Wavelet2_Default = Wavelet2<TSigma>;

static_assert(checkRankVector<Wavelet2_Default>);

}
