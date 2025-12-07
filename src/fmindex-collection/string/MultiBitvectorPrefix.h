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
struct MultiBitvectorPrefix {
    static constexpr size_t Sigma = TSigma;

    std::array<Bitvector, TSigma-1> bitvectors{};

    MultiBitvectorPrefix() = default;

    MultiBitvectorPrefix(std::span<uint8_t const> _symbols) {
        for (size_t sym{0}; sym < Sigma-1; ++sym) {
            bitvectors[sym] = Bitvector(std::views::iota(size_t{}, _symbols.size())
                                        | std::views::transform([&](size_t idx) {
                                            return _symbols[idx] <= sym; //!< or <=?
                                        }));
        }

        for (size_t i{0}; i < size(); ++i) {
            for (size_t sym{0}; sym < Sigma-2; ++sym) {
                auto v1 = bitvectors[sym].rank(i);
                auto v2 = bitvectors[sym+1].rank(i);
                (void)v1;
                (void)v2;
                assert(v1 <= v2);
            }
        }
    }

    size_t size() const {
        return bitvectors[0].size();
    }

    auto symbol(uint64_t idx) const -> uint64_t {
        assert(idx < size());
        for (size_t sym{0}; sym < Sigma-1; ++sym) {
            if (bitvectors[sym].symbol(idx)) {
                return sym;
            }
        }
        return Sigma-1;
    }

    uint64_t rank(uint64_t idx, uint64_t symb) const {
        assert(symb < TSigma);
        assert(idx <= size());
        auto v1 = prefix_rank(idx, symb+1);
        auto v2 = prefix_rank(idx, symb);
        auto v = v1 - v2;
        assert(v <= idx);
        return v;
    }

    uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
        assert(symb <= TSigma);
        assert(idx <= size());
        if (symb == 0) return 0;
        if (symb == Sigma) return idx;
        auto v = bitvectors[symb-1].rank(idx);
        assert(v <= idx);
        return v;
    }

    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
        assert(idx <= size());
        auto rs = std::array<uint64_t, TSigma>{};
        for (size_t sym{0}; sym < Sigma; ++sym) {
            rs[sym] = prefix_rank(idx, sym);
        }
        for (size_t sym{0}; sym < Sigma-1; ++sym) {
            rs[sym] = rs[sym+1] - rs[sym];
        }
        rs[Sigma-1] = idx - rs[Sigma-1];
        return rs;
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, TSigma>, std::array<uint64_t, TSigma>> {
        assert(idx <= size());

        auto rs = std::array<uint64_t, TSigma>{};
        auto prs = std::array<uint64_t, TSigma>{};
        for (size_t sym{0}; sym < Sigma; ++sym) {
            prs[sym] = prefix_rank(idx, sym);
        }
        for (size_t sym{0}; sym < Sigma-1; ++sym) {
            rs[sym] = prs[sym+1] - prs[sym];
        }
        rs[Sigma-1] = idx - prs[Sigma-1];
        return {rs, prs};
    }


    template <typename Archive>
    void serialize(this auto&& self, Archive& ar) {
        ar(self.bitvectors);
    }
};

template <uint64_t TSigma>
using MultiBitvectorPrefix_Bitvector = MultiBitvectorPrefix<TSigma>;

static_assert(checkString_c<MultiBitvectorPrefix_Bitvector>);

}
