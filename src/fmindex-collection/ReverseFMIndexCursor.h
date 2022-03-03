#pragma once

#include "ReverseFMIndex.h"

namespace fmindex_collection {

template <typename Index>
struct ReverseFMIndexCursor {
    static size_t constexpr Sigma = Index::Sigma;

    Index const* index;
    size_t lb;
    size_t len;

    ReverseFMIndexCursor()
        : index{nullptr}
    {}

    ReverseFMIndexCursor(Index const& index)
        : ReverseFMIndexCursor{index, 0, index.size()}
    {}

    ReverseFMIndexCursor(Index const& index, size_t lb, size_t len)
        : index{&index}
        , lb{lb}
        , len{len}
    {}
    ReverseFMIndexCursor(ReverseFMIndexCursor const&) = default;
    auto operator=(ReverseFMIndexCursor const&) -> ReverseFMIndexCursor& = default;

    auto extendLeft(uint8_t symb) const -> ReverseFMIndexCursor {
        size_t newLb  = index->occ.rank(lb, symb);
        size_t newLen = index->occ.rank(lb+len, symb) - newLb;
        return {*index, newLb, newLen};
    }
    auto extendLeft() const -> std::array<ReverseFMIndexCursor, Sigma> {
        auto [rs1, prs1] = index->occ.all_ranks(lb);
        auto [rs2, prs2] = index->occ.all_ranks(lb+len);

        auto cursors = std::array<ReverseFMIndexCursor, Sigma>{};
        cursors[0] = ReverseFMIndexCursor{*index, rs1[0], rs2[0] - rs1[0]};
        for (size_t i{1}; i < Sigma; ++i) {
            cursors[i] = ReverseFMIndexCursor{*index, rs1[i], rs2[i] - rs1[i]};
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
auto begin(ReverseFMIndexCursor<Index> const& _cursor) {
    return _cursor.lb;
}
template <typename Index>
auto end(ReverseFMIndexCursor<Index> const& _cursor) {
    return _cursor.lb + _cursor.len;
}


}
