// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "KMerFMIndex.h"

#include <compare>

namespace fmindex_collection {

template <typename Index>
struct KMerFMIndexCursor {
    static constexpr size_t Sigma    = Index::Sigma;
    static constexpr bool   Reversed = false;

    Index const* index{};
    size_t lb;
    size_t len{};

    KMerFMIndexCursor() noexcept = default;

    KMerFMIndexCursor(Index const& index) noexcept
        : KMerFMIndexCursor{index, 0, index.size()}
    {}

    KMerFMIndexCursor(Index const& index, size_t lb, size_t len) noexcept
        : index{&index}
        , lb{lb}
        , len{len}
    {}

    auto extendLeft(uint8_t symb) const -> KMerFMIndexCursor {
        size_t newLb  = index->bwt.rank(lb, symb);
        size_t newLen = index->bwt.rank(lb+len, symb) - newLb;
        return {*index, newLb + index->C[symb], newLen};
    }

    auto extendLeft() const -> std::array<KMerFMIndexCursor, Sigma> {
        auto [rs1, prs1] = index->bwt.all_ranks_and_prefix_ranks(lb);
        auto [rs2, prs2] = index->bwt.all_ranks_and_prefix_ranks(lb+len);

        for (size_t i{0}; i < rs1.size(); ++i) {
            rs1[i] += index->C[i];
            rs2[i] += index->C[i];
        }

        auto cursors = std::array<KMerFMIndexCursor, Sigma>{};
        cursors[0] = KMerFMIndexCursor{*index, rs1[0], rs2[0] - rs1[0]};
        for (size_t i{1}; i < Sigma; ++i) {
            cursors[i] = KMerFMIndexCursor{*index, rs1[i], rs2[i] - rs1[i]};
        }
        return cursors;
    }

    auto clipToKMer() const -> KMerFMIndexCursor {
        auto nlb = index->kmerStarts.gotoMarkingBwd(lb);
        auto nrb = index->kmerStarts.gotoMarkingFwd(lb+len);
        return KMerFMIndexCursor{*index, nlb, nrb-nlb};
    }

    bool empty() const {
        return len == 0;
    }

    size_t count() const {
        return len;
    }

};

template <typename Index>
auto begin(KMerFMIndexCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb};
}
template <typename Index>
auto end(KMerFMIndexCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb + _cursor.len};
}

}
