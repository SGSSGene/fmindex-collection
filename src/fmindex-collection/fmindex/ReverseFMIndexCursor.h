// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "ReverseFMIndex.h"

namespace fmindex_collection {

template <typename Index>
struct ReverseFMIndexCursor {
    static constexpr size_t Sigma    = Index::Sigma;
    static constexpr bool   Reversed = true;

    Index const* index{};
    size_t lb;
    size_t len{};
    size_t depth{};

    ReverseFMIndexCursor() noexcept = default;

    ReverseFMIndexCursor(Index const& index) noexcept
        : ReverseFMIndexCursor{index, 0, index.size(), 0}
    {}

    ReverseFMIndexCursor(Index const& index, size_t lb, size_t len, size_t depth) noexcept
        : index{&index}
        , lb{lb}
        , len{len}
        , depth{depth}
    {}

    auto extendRight(uint8_t symb) const -> ReverseFMIndexCursor {
        size_t newLb  = index->occ.rank(lb, symb);
        size_t newLen = index->occ.rank(lb+len, symb) - newLb;
        return {*index, newLb, newLen, depth+1};
    }
    auto extendRight() const -> std::array<ReverseFMIndexCursor, Sigma> {
        auto [rs1, prs1] = index->occ.all_ranks(lb);
        auto [rs2, prs2] = index->occ.all_ranks(lb+len);

        auto cursors = std::array<ReverseFMIndexCursor, Sigma>{};
        cursors[0] = ReverseFMIndexCursor{*index, rs1[0], rs2[0] - rs1[0], depth+1};
        for (size_t i{1}; i < Sigma; ++i) {
            cursors[i] = ReverseFMIndexCursor{*index, rs1[i], rs2[i] - rs1[i], depth+1};
        }
        return cursors;
    }

    bool empty() const {
        return len == 0;
    }

    size_t count() const {
        return len;
    }

    size_t query_length() const {
        return depth;
    }

};

template <typename Index>
auto begin(ReverseFMIndexCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb};
}
template <typename Index>
auto end(ReverseFMIndexCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb + _cursor.len};
}


}
