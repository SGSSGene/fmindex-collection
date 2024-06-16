// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "concepts.h"

#include <vector>

namespace fmindex_collection::rankvector {

template <size_t TSigma>
struct Naive {
    static constexpr size_t Sigma = TSigma;

    std::vector<std::array<uint64_t, Sigma>> occ{};
    size_t                                   totalLength;

    Naive() = default;
    Naive(std::span<uint8_t const> _symbols) {
        occ.reserve(_symbols.size()+1);
        occ.push_back(std::array<uint64_t, Sigma>{}); // Initial first row with only zeros

        for (uint64_t i{0}; i < _symbols.size(); ++i) {
            occ.push_back(occ.back());
            occ.back()[_symbols[i]] += 1;
        }
        totalLength = _symbols.size();
    }

    size_t size() const noexcept {
        return totalLength;
    }

    uint8_t symbol(size_t idx) const noexcept {
        idx += 1;
        for (uint64_t i{0}; i < Sigma-1; ++i) {
            if (occ[idx][i] > occ[idx-1][i]) {
                return i;
            }
        }
        return Sigma-1;
    }

    uint64_t rank(size_t idx, uint8_t symb) const noexcept {
        return occ[idx][symb];
    }

    uint64_t prefix_rank(size_t idx, uint8_t symb) const noexcept {
        uint64_t a{};
        for (uint64_t i{0}; i <= symb; ++i) {
            a += occ[idx][i];
        }
        return a;
    }

    auto all_ranks(size_t idx) const -> std::array<uint64_t, Sigma> {
        return occ[idx];
    }

    auto all_ranks_and_prefix_ranks(size_t idx) const -> std::tuple<std::array<uint64_t, Sigma>, std::array<uint64_t, Sigma>> {
        auto rs  = occ[idx];
        auto prs = rs;
        for (size_t i{1}; i < prs.size(); ++i) {
            prs[i] = prs[i] + prs[i-1];
        }
        return {rs, prs};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(occ, totalLength);
    }
};
static_assert(checkRankVector<Naive>);

}
