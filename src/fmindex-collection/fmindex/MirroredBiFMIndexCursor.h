// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "MirroredBiFMIndex.h"

namespace fmindex_collection {

template <typename Index>
struct LeftMirroredBiFMIndexCursor;

template <typename Index>
struct MirroredBiFMIndexCursor {
    static constexpr size_t Sigma    = Index::Sigma;
    static constexpr bool   Reversed = false;

    Index const* index{};
    size_t lb;
    size_t lbRev;
    size_t len{};
    MirroredBiFMIndexCursor() noexcept = default;
    MirroredBiFMIndexCursor(Index const& index) noexcept
        : MirroredBiFMIndexCursor{index, 0, 0, index.size()}
    {}
    MirroredBiFMIndexCursor(Index const& index, size_t lb, size_t lbRev, size_t len) noexcept
        : index{&index}
        , lb{lb}
        , lbRev{lbRev}
        , len{len}
    {}

    bool operator==(MirroredBiFMIndexCursor const& _other) const noexcept {
        return lb == _other.lb
               && len == _other.len;
    }
    bool empty() const {
        return len == 0;
    }
    size_t count() const {
        return len;
    }
    auto extendLeft() const -> std::array<MirroredBiFMIndexCursor, Sigma> {
        auto const& bwt = index->bwt;
        auto [rs1, prs1] = bwt.all_ranks_and_prefix_ranks(lb);
        auto [rs2, prs2] = bwt.all_ranks_and_prefix_ranks(lb+len);

        for (size_t i{0}; i < rs1.size(); ++i) {
            rs1[i] += index->C[i];
            rs2[i] += index->C[i];
        }

        auto cursors = std::array<MirroredBiFMIndexCursor, Sigma>{};
        for (size_t i{0}; i < Sigma; ++i) {
            cursors[i] = MirroredBiFMIndexCursor{*index, rs1[i], lbRev + prs2[i] - prs1[i], rs2[i] - rs1[i]};
        }
        return cursors;
    }

    auto extendRight() const -> std::array<MirroredBiFMIndexCursor, Sigma> {
        auto const& bwt = index->bwt;
        auto [rs1, prs1] = bwt.all_ranks_and_prefix_ranks(lbRev);
        auto [rs2, prs2] = bwt.all_ranks_and_prefix_ranks(lbRev+len);

        for (size_t i{0}; i < rs1.size(); ++i) {
            rs1[i] += index->C[i];
            rs2[i] += index->C[i];
        }

        auto cursors = std::array<MirroredBiFMIndexCursor, Sigma>{};
        for (size_t i{0}; i < Sigma; ++i) {
            cursors[i] = MirroredBiFMIndexCursor{*index, lb + prs2[i] - prs1[i], rs1[i], rs2[i] - rs1[i]};
        }
        return cursors;
    }
    void prefetchLeft() const {
/*        if constexpr (OccTablePrefetch<Index>) {
            auto& occ = index->occ;
            occ.prefetch(lb);
            occ.prefetch(lb+len);
        }*/
    }
    void prefetchRight() const {
/*        if constexpr (OccTablePrefetch<Index>) {
            auto& occ = index->occ;
            occ.prefetch(lbRev);
            occ.prefetch(lbRev+len);
        }*/
    }

    auto extendLeft(size_t symb) const -> MirroredBiFMIndexCursor {
        auto& bwt = index->bwt;
        size_t newLb    = bwt.rank(lb, symb);
        size_t newLbRev = lbRev + bwt.prefix_rank(lb+len, symb) - bwt.prefix_rank(lb, symb);
        size_t newLen   = bwt.rank(lb+len, symb) - newLb;
        auto newCursor = MirroredBiFMIndexCursor{*index, newLb + index->C[symb], newLbRev, newLen};
        return newCursor;
    }
    auto extendRight(size_t symb) const -> MirroredBiFMIndexCursor {
        auto& bwt = index->bwt;
        size_t newLb    = lb + bwt.prefix_rank(lbRev+len, symb) - bwt.prefix_rank(lbRev, symb);
        size_t newLbRev = bwt.rank(lbRev, symb);
        size_t newLen   = bwt.rank(lbRev+len, symb) - newLbRev;
        auto newCursor = MirroredBiFMIndexCursor{*index, newLb, newLbRev + index->C[symb], newLen};
        return newCursor;
    }
};

template <typename Index>
auto begin(MirroredBiFMIndexCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb};
}
template <typename Index>
auto end(MirroredBiFMIndexCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb + _cursor.len};
}

template <typename Index>
struct LeftMirroredBiFMIndexCursor {
    static constexpr size_t Sigma    = Index::Sigma;
    static constexpr bool   Reversed = false;

    Index const* index;
    size_t lb;
    size_t len;
    LeftMirroredBiFMIndexCursor(MirroredBiFMIndexCursor<Index> const& _other)
        : index{_other.index}
        , lb{_other.lb}
        , len{_other.len}
    {}
    LeftMirroredBiFMIndexCursor()
        : index{nullptr}
    {}
    LeftMirroredBiFMIndexCursor(Index const& index)
        : LeftMirroredBiFMIndexCursor{index, 0, index.size()}
    {}
    LeftMirroredBiFMIndexCursor(Index const& index, size_t lb, size_t len)
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
    auto extendLeft() const -> std::array<LeftMirroredBiFMIndexCursor, Sigma> {
        auto const& bwt = index->bwt;
        auto [rs1, prs1] = bwt.all_ranks_and_prefix_ranks(lb);
        auto [rs2, prs2] = bwt.all_ranks_and_prefix_ranks(lb+len);

        for (size_t i{0}; i < rs1.size(); ++i) {
            rs1[i] += index->C[i];
            rs2[i] += index->C[i];
        }

        auto cursors = std::array<LeftMirroredBiFMIndexCursor, Sigma>{};
        cursors[0] = LeftMirroredBiFMIndexCursor{*index, rs1[0], rs2[0] - rs1[0]};
        for (size_t i{1}; i < Sigma; ++i) {
            cursors[i] = LeftMirroredBiFMIndexCursor{*index, rs1[i], rs2[i] - rs1[i]};
        }
        return cursors;
    }

    auto extendLeft(size_t symb) const -> LeftMirroredBiFMIndexCursor {
        auto& bwt = index->bwt;

        size_t newLb    = bwt.rank(lb, symb);
        size_t newLen   = bwt.rank(lb+len, symb) - newLb;

        auto newCursor = LeftMirroredBiFMIndexCursor{*index, newLb + index->C[symb], newLen};
        return newCursor;
    }
};

template <typename Index>
auto begin(LeftMirroredBiFMIndexCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb};
}
template <typename Index>
auto end(LeftMirroredBiFMIndexCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb + _cursor.len};
}

}

namespace std {

template <typename index_t>
struct hash<fmindex_collection::MirroredBiFMIndexCursor<index_t>> {
    auto operator()(fmindex_collection::MirroredBiFMIndexCursor<index_t> const& cursor) const -> size_t {
        return hash<size_t>()(cursor.lb)
            ^ hash<size_t>()(cursor.len);
    }
};

}
