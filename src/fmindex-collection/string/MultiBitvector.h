// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../bitvector/Bitvector.h"
#include "concepts.h"
#include "utils.h"

#include <bitset>
#include <ranges>
#include <vector>

namespace fmc::string {

template <size_t TSigma, Bitvector_c Bitvector = ::fmc::bitvector::Bitvector>
struct MultiBitvector {
    static constexpr size_t Sigma = TSigma;

    std::array<Bitvector, TSigma> bitvectors{};

    MultiBitvector() = default;

    MultiBitvector(std::span<uint8_t const> _symbols) {
        for (size_t sym{0}; sym < Sigma; ++sym) {
            bitvectors[sym] = Bitvector(std::views::iota(size_t{}, _symbols.size())
                                        | std::views::transform([&](size_t idx) {
                                            return _symbols[idx] == sym;
                                        }));
        }
    }

    void prefetch(uint64_t idx) const {
        if constexpr (requires() { {bitvectors[0].prefetch(idx)}; }) {
            for (auto const& bv : bitvectors) {
                bv.prefetch(idx);
            }
        }
    }

    size_t size() const {
        return bitvectors[0].size();
    }

    uint8_t symbol(uint64_t idx) const {
        assert(idx < size());
        for (size_t sym{0}; sym < Sigma; ++sym) {
            if (bitvectors[sym].symbol(idx)) {
                return sym;
            }
        }
        return 0;
    }

    uint64_t rank(uint64_t idx, uint8_t symb) const {
        assert(symb < TSigma);
        assert(idx <= size());

        return bitvectors[symb].rank(idx);
    }

    uint64_t prefix_rank(uint64_t idx, uint8_t symb) const {
        assert(symb <= TSigma);
        assert(idx <= size());
        size_t a{};
        for (size_t i{1}; i <= symb; ++i) {
            a += bitvectors[i-1].rank(idx);
        }
        return a;
    }

    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
        assert(idx <= size());
        auto rs = std::array<uint64_t, TSigma>{};
        for (size_t sym{0}; sym < Sigma; ++sym) {
            rs[sym] = rank(idx, sym);
        }
        return rs;
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, TSigma>, std::array<uint64_t, TSigma>> {
        assert(idx <= size());

        auto rs  = all_ranks(idx);
        auto prs = std::array<uint64_t, TSigma>{};
        for (size_t i{1}; i < Sigma; ++i) {
            prs[i] = prs[i-1] + rs[i-1];
        }
        return {rs, prs};
    }


    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bitvectors);
    }
};

template <uint64_t TSigma>
using MultiBitvector_Bitvector = MultiBitvector<TSigma>;

static_assert(checkString_c<MultiBitvector_Bitvector>);

}
