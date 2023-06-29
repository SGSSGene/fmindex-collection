// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "BiFMIndex.h"

namespace fmindex_collection {

template <typename Index>
struct LeftBiFMIndexCursor;

template <typename Index>
struct BiFMIndexCursor {
    static constexpr size_t Sigma    = Index::Sigma;
    static constexpr bool   Reversed = false;

    Index const* index{};
    size_t lb;
    size_t lbRev;
    size_t len{};
    BiFMIndexCursor() noexcept = default;
    BiFMIndexCursor(Index const& index) noexcept
        : BiFMIndexCursor{index, 0, 0, index.size()}
    {}
    BiFMIndexCursor(Index const& index, size_t lb, size_t lbRev, size_t len) noexcept
        : index{&index}
        , lb{lb}
        , lbRev{lbRev}
        , len{len}
    {}

    bool operator==(BiFMIndexCursor const& _other) const noexcept {
        return lb == _other.lb
               && len == _other.len;
    }
    bool empty() const {
        return len == 0;
    }
    size_t count() const {
        return len;
    }
    auto extendLeft() const -> std::array<BiFMIndexCursor, Sigma> {
        auto const& occ = index->occ;
        if constexpr (OccTablePrefetch<Index>) {
            occ.prefetch(lb+len);
        }
        auto [rs1, prs1] = occ.all_ranks(lb);
        auto [rs2, prs2] = occ.all_ranks(lb+len);

        auto cursors = std::array<BiFMIndexCursor, Sigma>{};
        cursors[0] = BiFMIndexCursor{*index, rs1[0], lbRev, rs2[0] - rs1[0]};
        cursors[0].prefetchLeft();
        for (size_t i{1}; i < Sigma; ++i) {
            cursors[i] = BiFMIndexCursor{*index, rs1[i], lbRev + prs2[i-1] - prs1[i-1], rs2[i] - rs1[i]};
        }
        return cursors;
    }

    auto extendRight() const -> std::array<BiFMIndexCursor, Sigma> {
        auto const& occ = index->occRev;
        if constexpr (OccTablePrefetch<Index>) {
            occ.prefetch(lbRev+len);
        }
        auto [rs1, prs1] = occ.all_ranks(lbRev);
        auto [rs2, prs2] = occ.all_ranks(lbRev+len);

        auto cursors = std::array<BiFMIndexCursor, Sigma>{};
        cursors[0] = BiFMIndexCursor{*index, lb, rs1[0], rs2[0] - rs1[0]};
        cursors[0].prefetchRight();
        for (size_t i{1}; i < Sigma; ++i) {
            cursors[i] = BiFMIndexCursor{*index, lb + prs2[i-1] - prs1[i-1], rs1[i], rs2[i] - rs1[i]};
        }
        return cursors;
    }
    void prefetchLeft() const {
        if constexpr (OccTablePrefetch<Index>) {
            auto& occ = index->occ;
            occ.prefetch(lb);
            occ.prefetch(lb+len);
        }
    }
    void prefetchRight() const {
        if constexpr (OccTablePrefetch<Index>) {
            auto& occ = index->occRev;
            occ.prefetch(lbRev);
            occ.prefetch(lbRev+len);
        }
    }

    auto extendLeft(size_t symb) const -> BiFMIndexCursor {
        assert(symb > 0);
        auto& occ = index->occ;
        size_t newLb    = occ.rank(lb, symb);
        size_t newLbRev = lbRev + occ.prefix_rank(lb+len, symb-1) - occ.prefix_rank(lb, symb-1);
        size_t newLen   = occ.rank(lb+len, symb) - newLb;
        auto newCursor = BiFMIndexCursor{*index, newLb, newLbRev, newLen};
        newCursor.prefetchLeft();
        return newCursor;
    }
    auto extendRight(size_t symb) const -> BiFMIndexCursor {
        assert(symb > 0);
        auto& occ = index->occRev;
        size_t newLb    = lb + occ.prefix_rank(lbRev+len, symb-1) - occ.prefix_rank(lbRev, symb-1);
        size_t newLbRev = occ.rank(lbRev, symb);
        size_t newLen   = occ.rank(lbRev+len, symb) - newLbRev;
        auto newCursor = BiFMIndexCursor{*index, newLb, newLbRev, newLen};
        newCursor.prefetchRight();
        return newCursor;
    }
};

template <typename Index>
auto begin(BiFMIndexCursor<Index> const& _cursor) {
    return _cursor.lb;
}
template <typename Index>
auto end(BiFMIndexCursor<Index> const& _cursor) {
    return _cursor.lb + _cursor.len;
}

template <typename Index>
struct LeftBiFMIndexCursor {
    static constexpr size_t Sigma    = Index::Sigma;
    static constexpr bool   Reversed = false;

    Index const* index;
    size_t lb;
    size_t len;
    LeftBiFMIndexCursor(BiFMIndexCursor<Index> const& _other)
        : index{_other.index}
        , lb{_other.lb}
        , len{_other.len}
    {}
    LeftBiFMIndexCursor()
        : index{nullptr}
    {}
    LeftBiFMIndexCursor(Index const& index)
        : LeftBiFMIndexCursor{index, 0, index.size()}
    {}
    LeftBiFMIndexCursor(Index const& index, size_t lb, size_t len)
        : index{&index}
        , lb{lb}
        , len{len}
    {}
    bool empty() const {
        return len == 0;
    }
    size_t count() const {
        return len;
    }
    auto extendLeft() const -> std::array<LeftBiFMIndexCursor, Sigma> {
        auto const& occ = index->occ;
        auto [rs1, prs1] = occ.all_ranks(lb);
        auto [rs2, prs2] = occ.all_ranks(lb+len);

        auto cursors = std::array<LeftBiFMIndexCursor, Sigma>{};
        cursors[0] = LeftBiFMIndexCursor{*index, rs1[0], rs2[0] - rs1[0]};
        for (size_t i{1}; i < Sigma; ++i) {
            cursors[i] = LeftBiFMIndexCursor{*index, rs1[i], rs2[i] - rs1[i]};
//            cursors[i].prefetchLeft();
        }
        return cursors;
    }

    auto extendLeft(size_t symb) const -> LeftBiFMIndexCursor {
        assert(symb > 0);
        auto& occ = index->occ;

        size_t newLb    = occ.rank(lb, symb);
        size_t newLen   = occ.rank(lb+len, symb) - newLb;
        if constexpr (OccTablePrefetch<Index>) {
            occ.prefetch(newLb);
            occ.prefetch(newLb + newLen);
        }

        auto newCursor = LeftBiFMIndexCursor{*index, newLb, newLen};
        return newCursor;
    }
};
}

namespace std {

template <typename index_t>
struct hash<fmindex_collection::BiFMIndexCursor<index_t>> {
    auto operator()(fmindex_collection::BiFMIndexCursor<index_t> const& cursor) const -> size_t {
        return hash<size_t>()(cursor.lb)
            ^ hash<size_t>()(cursor.len);
    }
};

}
