// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "LinearFMIndex.h"

#include <compare>

namespace fmc {

template <typename Index>
struct LinearFMIndexCursor {
    static constexpr size_t Sigma    = Index::Sigma;
    static constexpr bool   Reversed = false;

    Index const* index{};
    size_t lb;
    size_t len{};
    size_t col{};

    LinearFMIndexCursor() noexcept = default;

    LinearFMIndexCursor(Index const& index) noexcept
        : LinearFMIndexCursor{index, 0, index.size(), index.columns.size()-1}
    {}

    LinearFMIndexCursor(Index const& index, size_t lb, size_t len, size_t col) noexcept
        : index{&index}
        , lb{lb}
        , len{len}
        , col{col}
    {}

    auto extendLeft(uint8_t symb) const -> LinearFMIndexCursor {
        auto const& column = index->columns[col];
        size_t newLb  = column.bwt.rank(lb, symb);
        size_t newLen = column.bwt.rank(lb+len, symb) - newLb;
        return {*index, newLb + column.C[symb], newLen, col-1};
    }
    auto extendLeft() const -> std::array<LinearFMIndexCursor, Sigma> {
        auto const& column = index->columns[col];
        auto [rs1, prs1] = column.bwt.all_ranks_and_prefix_ranks(lb);
        auto [rs2, prs2] = column.bwt.all_ranks_and_prefix_ranks(lb+len);

        for (size_t i{0}; i < rs1.size(); ++i) {
            rs1[i] += column.C[i];
            rs2[i] += column.C[i];
        }

        auto cursors = std::array<LinearFMIndexCursor, Sigma>{};
        cursors[0] = LinearFMIndexCursor{*index, rs1[0], rs2[0] - rs1[0], col-1};
        for (size_t i{1}; i < Sigma; ++i) {
            cursors[i] = LinearFMIndexCursor{*index, rs1[i], rs2[i] - rs1[i], col-1};
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
auto begin(LinearFMIndexCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb};
}
template <typename Index>
auto end(LinearFMIndexCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb + _cursor.len};
}

}
