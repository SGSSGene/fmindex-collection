// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "FMIndex.h"

namespace fmindex_collection {

template <typename Index>
struct FMIndexCursor {
    static constexpr size_t Sigma    = Index::Sigma;
    static constexpr bool   Reversed = false;

    Index const* index{};
    size_t lb;
    size_t len{};

    FMIndexCursor() noexcept = default;

    FMIndexCursor(Index const& index) noexcept
        : FMIndexCursor{index, 0, index.size()}
    {}

    FMIndexCursor(Index const& index, size_t lb, size_t len) noexcept
        : index{&index}
        , lb{lb}
        , len{len}
    {}

    auto extendLeft(uint8_t symb) const -> FMIndexCursor {
        size_t newLb  = index->occ.rank(lb, symb);
        size_t newLen = index->occ.rank(lb+len, symb) - newLb;
        return {*index, newLb, newLen};
    }
    auto extendLeft() const -> std::array<FMIndexCursor, Sigma> {
        auto [rs1, prs1] = index->occ.all_ranks(lb);
        auto [rs2, prs2] = index->occ.all_ranks(lb+len);

        auto cursors = std::array<FMIndexCursor, Sigma>{};
        cursors[0] = FMIndexCursor{*index, rs1[0], rs2[0] - rs1[0]};
        for (size_t i{1}; i < Sigma; ++i) {
            cursors[i] = FMIndexCursor{*index, rs1[i], rs2[i] - rs1[i]};
        }
        return cursors;
    }

    bool empty() const {
        return len == 0;
    }

    size_t count() const {
        return len;
    }

};

template <typename Index>
auto begin(FMIndexCursor<Index> const& _cursor) {
    return _cursor.lb;
}
template <typename Index>
auto end(FMIndexCursor<Index> const& _cursor) {
    return _cursor.lb + _cursor.len;
}


}
