// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "BiFMIndex.h"

namespace fmc {

template <typename Index>
struct LeftBiFMIndexCursor;

template <typename Index>
struct BiFMIndexCursor {
    static constexpr size_t Sigma    = Index::Sigma;
    static constexpr bool   Reversed = false;

    constexpr bool static HasDualRank = requires(Index::String str, size_t idx) {
        { str.all_ranks_dual(idx, idx, [](size_t, size_t, size_t, size_t, size_t) {}) };
    };

    Index const* index{};
    size_t lb;
    size_t lbRev;
    size_t len{};
    size_t steps{}; // number of extension steps taken
    BiFMIndexCursor() noexcept = default;
    BiFMIndexCursor(Index const& index) noexcept
        : BiFMIndexCursor{index, 0, 0, index.size(), 0}
    {}
    BiFMIndexCursor(Index const& index, size_t lb, size_t lbRev, size_t len, size_t steps) noexcept
        : index{&index}
        , lb{lb}
        , lbRev{lbRev}
        , len{len}
        , steps{steps}
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

    auto fetchRightBwt() const -> auto const& {
        if constexpr (std::same_as<typename Index::RevBwtType, std::nullptr_t>) {
            return index->bwt;
        } else {
            return index->bwtRev;
        }
    }

    auto extendLeft() const -> std::array<BiFMIndexCursor, Sigma> requires (!HasDualRank) {
        auto cursors = std::array<BiFMIndexCursor, Sigma>{};

        auto const& bwt = index->bwt;
        auto [rs1, prs1] = bwt.all_ranks_and_prefix_ranks(lb);
        auto [rs2, prs2] = bwt.all_ranks_and_prefix_ranks(lb+len);

        for (size_t i{0}; i < Sigma; ++i) {
            cursors[i] = BiFMIndexCursor{*index, rs1[i] + index->C[i], lbRev + prs2[i] - prs1[i], rs2[i] - rs1[i], steps+1};
        }
        return cursors;
    }

    auto extendRight() const -> std::array<BiFMIndexCursor, Sigma> requires (!HasDualRank) {
        auto cursors = std::array<BiFMIndexCursor, Sigma>{};

        auto const& bwt = fetchRightBwt();
        auto [rs1, prs1] = bwt.all_ranks_and_prefix_ranks(lbRev);
        auto [rs2, prs2] = bwt.all_ranks_and_prefix_ranks(lbRev+len);

        for (size_t i{0}; i < Sigma; ++i) {
            cursors[i] = BiFMIndexCursor{*index, lb + prs2[i] - prs1[i], rs1[i] + index->C[i], rs2[i] - rs1[i], steps+1};
        }
        return cursors;
    }

    auto extendLeft() const -> std::array<BiFMIndexCursor, Sigma> requires HasDualRank {
        auto ret = std::array<BiFMIndexCursor, Sigma>{};
        auto& bwt = index->bwt;
        bwt.all_ranks_dual(lb, lb+len, [&](size_t symb, size_t rs1, size_t rs2, size_t prs1, size_t prs2) {
            auto newLb    = index->C[symb] + rs1;
            auto newLen   = rs2 - rs1;
            auto newLbRev = lbRev + prs2 - prs1;
            ret[symb] = BiFMIndexCursor{*index, newLb, newLbRev, newLen, steps+1};
        });
        return ret;
    }
    auto extendRight() const -> std::array<BiFMIndexCursor, Sigma> requires HasDualRank {
        auto ret = std::array<BiFMIndexCursor, Sigma>{};
        auto& bwt = fetchRightBwt();
        bwt.all_ranks_dual(lbRev, lbRev+len, [&](size_t symb, size_t rs1, size_t rs2, size_t prs1, size_t prs2) {
            auto newLbRev = index->C[symb] + rs1;
            auto newLen   = rs2 - rs1;
            auto newLb    = lb + prs2 - prs1;
            ret[symb] = BiFMIndexCursor{*index, newLb, newLbRev, newLen, steps+1};
        });
        return ret;
    }


    void prefetchLeft() const {
    }
    void prefetchRight() const {
    }

    auto extendLeft(size_t symb) const -> BiFMIndexCursor {
        auto& bwt = index->bwt;
        size_t newLb    = bwt.rank(lb, symb);
        size_t newLbRev = lbRev + bwt.prefix_rank(lb+len, symb) - bwt.prefix_rank(lb, symb);
        size_t newLen   = bwt.rank(lb+len, symb) - newLb;
        auto newCursor = BiFMIndexCursor{*index, newLb + index->C[symb], newLbRev, newLen, steps+1};
        return newCursor;
    }
    auto extendRight(size_t symb) const -> BiFMIndexCursor {
        auto const& bwt = fetchRightBwt();
        size_t newLb    = lb + bwt.prefix_rank(lbRev+len, symb) - bwt.prefix_rank(lbRev, symb);
        size_t newLbRev = bwt.rank(lbRev, symb);
        size_t newLen   = bwt.rank(lbRev+len, symb) - newLbRev;
        auto newCursor = BiFMIndexCursor{*index, newLb, newLbRev + index->C[symb], newLen, steps+1};
        return newCursor;
    }

    // This requires that all rows have the same BWT entry (or only a single one is available)
    // - must have at least marked a single row
    // - all rows must have the same 'bwt' symbol
    auto extendLeftBySymbol(size_t symb) const -> BiFMIndexCursor {
        auto& bwt = index->bwt;

        assert(count() > 0);
        assert([&]() {
            for (size_t i{lb}; i < lb + len; ++i) {
                if (symb != bwt.symbol(i)) {
                    return false;
                }
            }
            return true;
        }());

        size_t newLb    = bwt.rank(lb, symb);
        size_t newLbRev = lbRev;
        size_t newLen   = len;
        auto newCursor  = BiFMIndexCursor{*index, newLb + index->C[symb], newLbRev, newLen, steps+1};
        return newCursor;
    }

    // see extendLeftBySymbol
    auto extendRightBySymbol(size_t symb) const -> BiFMIndexCursor {
        auto const& bwt = fetchRightBwt();
        size_t newLb    = lb;
        size_t newLbRev = bwt.rank(lbRev, symb);
        size_t newLen   = len;
        auto newCursor = BiFMIndexCursor{*index, newLb, newLbRev + index->C[symb], newLen, steps+1};
        return newCursor;
    }

    // This requires that all rows have the same BWT entry (or only a single one is available)
    // - must have at least marked a single row
    // - all rows must have the same 'bwt' symbol
    auto extendLeftBySymbol() const -> std::tuple<size_t, BiFMIndexCursor> {
        auto& bwt = index->bwt;
        auto symb = bwt.symbol(lb);
        return {symb, extendLeftBySymbol(symb)};
    }

    // see extendLeftBySymbol
    auto extendRightBySymbol() const -> std::tuple<size_t, BiFMIndexCursor> {
        auto const& bwt = fetchRightBwt();
        auto symb = bwt.symbol(lbRev);
        return {symb, extendRightBySymbol(symb)};
    }


    auto symbolLeft() const -> size_t {
        return index->bwt.symbol(lb);
    }
    auto symbolRight() const -> size_t {
        // Reuse bwt if bwtrev is not available
        if constexpr (std::same_as<typename Index::RevBwtType, std::nullptr_t>) {
            return index->bwt.symbol(lbRev);
        } else {
            return index->bwtRev.symbol(lbRev);
        }
    }
};

template <typename Index>
auto begin(BiFMIndexCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb};
}
template <typename Index>
auto end(BiFMIndexCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb + _cursor.len};
}

template <typename Index>
struct LeftBiFMIndexCursor {
    static constexpr size_t Sigma    = Index::Sigma;
    static constexpr bool   Reversed = false;

    Index const* index;
    size_t lb;
    size_t len;
    size_t steps;
    LeftBiFMIndexCursor(BiFMIndexCursor<Index> const& _other)
        : index{_other.index}
        , lb{_other.lb}
        , len{_other.len}
        , steps{_other.steps}
    {}
    LeftBiFMIndexCursor()
        : index{nullptr}
    {}
    LeftBiFMIndexCursor(Index const& index)
        : LeftBiFMIndexCursor{index, 0, index.size(), 0}
    {}
    LeftBiFMIndexCursor(Index const& index, size_t lb, size_t len, size_t steps)
        : index{&index}
        , lb{lb}
        , len{len}
        , steps{steps}
    {}
    bool empty() const {
        return len == 0;
    }
    size_t count() const {
        return len;
    }
    auto extendLeft() const -> std::array<LeftBiFMIndexCursor, Sigma> {
        auto cursors = std::array<LeftBiFMIndexCursor, Sigma>{};

        auto const& bwt = index->bwt;
        auto [rs1, prs1] = bwt.all_ranks_and_prefix_ranks(lb);
        auto [rs2, prs2] = bwt.all_ranks_and_prefix_ranks(lb+len);

        for (size_t i{0}; i < Sigma; ++i) {
            cursors[i] = LeftBiFMIndexCursor{*index, rs1[i] + index->C[i], rs2[i] - rs1[i], steps+1};
        }
        return cursors;
    }

    auto extendLeft(size_t symb) const -> LeftBiFMIndexCursor {
        auto& bwt = index->bwt;

        size_t newLb    = bwt.rank(lb, symb);
        size_t newLen   = bwt.rank(lb+len, symb) - newLb;
        auto newCursor = LeftBiFMIndexCursor{*index, newLb + index->C[symb], newLen, steps+1};
        return newCursor;
    }
};

template <typename Index>
auto begin(LeftBiFMIndexCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb};
}
template <typename Index>
auto end(LeftBiFMIndexCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb + _cursor.len};
}

}

namespace std {

template <typename index_t>
struct hash<fmc::BiFMIndexCursor<index_t>> {
    auto operator()(fmc::BiFMIndexCursor<index_t> const& cursor) const -> size_t {
        return hash<size_t>()(cursor.lb)
            ^ hash<size_t>()(cursor.len);
    }
};

}
