// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "concepts.h"

#include "../bitvector/CompactBitvector.h"

#include <bitset>
#include <cassert>
#include <vector>

namespace fmindex_collection::rankvector {

struct CompactBitvector {
    static constexpr size_t Sigma = 2;

    bitvector::CompactBitvector bitvector;

    CompactBitvector(std::span<uint8_t const> _symbols)
        : CompactBitvector(_symbols, 1)
    {}

    CompactBitvector(std::span<uint8_t const> _symbols, uint8_t _symbol)
        : bitvector{_symbols.size(), [&](size_t idx) {
            return _symbols[idx] == _symbol;
        }}
    {}

    CompactBitvector() = default;
    CompactBitvector(CompactBitvector const&) = default;
    CompactBitvector(CompactBitvector&&) noexcept = default;
    auto operator=(CompactBitvector const&) -> CompactBitvector& = default;
    auto operator=(CompactBitvector&&) noexcept -> CompactBitvector& = default;


    size_t size() const noexcept {
        return bitvector.size();
    }

    uint8_t symbol(size_t idx) const noexcept {
        return bitvector.symbol(idx);
    }

    uint64_t rank(size_t idx, uint8_t symb=1) const noexcept {
        auto v = bitvector.rank(idx);
        if (symb == 0) {
            v = idx - v;
        }
        return v;
    }

    uint64_t prefix_rank(size_t idx, uint8_t symb) const noexcept {
        if (symb == 1) return idx;
        return rank(idx, 0);
    }

    auto all_ranks(size_t idx) const -> std::array<uint64_t, Sigma> {
        auto r = rank(idx);
        return {idx-r, r};
    }

    auto all_ranks_and_prefix_ranks(size_t idx) const -> std::tuple<std::array<uint64_t, Sigma>, std::array<uint64_t, Sigma>> {
        auto rs  = all_ranks(idx);
        auto prs = rs;
        prs[1] = idx;
        return {rs, prs};
    }


    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bitvector);
    }
};

static_assert(RankVector<CompactBitvector>);

}
