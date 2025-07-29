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
        auto const& bwt = index->bwt;
        auto [rs1, prs1] = bwt.all_ranks_and_prefix_ranks(lb);
        auto [rs2, prs2] = bwt.all_ranks_and_prefix_ranks(lb+len);

        for (size_t i{0}; i < rs1.size(); ++i) {
            rs1[i] += index->C[i];
            rs2[i] += index->C[i];
        }

        auto cursors = std::array<BiFMIndexCursor, Sigma>{};
        for (size_t i{0}; i < Sigma; ++i) {
            cursors[i] = BiFMIndexCursor{*index, rs1[i], lbRev + prs2[i] - prs1[i], rs2[i] - rs1[i]};
        }
        return cursors;
    }

    auto extendRight() const -> std::array<BiFMIndexCursor, Sigma> {
        // Reuse bwt if bwtrev is not available
        if constexpr (std::same_as<typename Index::RevBwtType, std::nullptr_t>) {
            auto const& bwt = index->bwt;
            auto [rs1, prs1] = bwt.all_ranks_and_prefix_ranks(lbRev);
            auto [rs2, prs2] = bwt.all_ranks_and_prefix_ranks(lbRev+len);

            for (size_t i{0}; i < rs1.size(); ++i) {
                rs1[i] += index->C[i];
                rs2[i] += index->C[i];
            }

            auto cursors = std::array<BiFMIndexCursor, Sigma>{};
            for (size_t i{0}; i < Sigma; ++i) {
                cursors[i] = BiFMIndexCursor{*index, lb + prs2[i] - prs1[i], rs1[i], rs2[i] - rs1[i]};
            }
            return cursors;
        } else {
            auto const& bwt = index->bwtRev;
            auto [rs1, prs1] = bwt.all_ranks_and_prefix_ranks(lbRev);
            auto [rs2, prs2] = bwt.all_ranks_and_prefix_ranks(lbRev+len);

            for (size_t i{0}; i < rs1.size(); ++i) {
                rs1[i] += index->C[i];
                rs2[i] += index->C[i];
            }

            auto cursors = std::array<BiFMIndexCursor, Sigma>{};
            for (size_t i{0}; i < Sigma; ++i) {
                cursors[i] = BiFMIndexCursor{*index, lb + prs2[i] - prs1[i], rs1[i], rs2[i] - rs1[i]};
            }
            return cursors;
        }
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
        auto newCursor = BiFMIndexCursor{*index, newLb + index->C[symb], newLbRev, newLen};
        return newCursor;
    }
    auto extendRight(size_t symb) const -> BiFMIndexCursor {
        if constexpr (std::same_as<typename Index::RevBwtType, std::nullptr_t>) {
            auto const& bwt = index->bwt;
            size_t newLb    = lb + bwt.prefix_rank(lbRev+len, symb) - bwt.prefix_rank(lbRev, symb);
            size_t newLbRev = bwt.rank(lbRev, symb);
            size_t newLen   = bwt.rank(lbRev+len, symb) - newLbRev;
            auto newCursor = BiFMIndexCursor{*index, newLb, newLbRev + index->C[symb], newLen};
            return newCursor;
        } else {
            auto const& bwt = index->bwtRev;
            size_t newLb    = lb + bwt.prefix_rank(lbRev+len, symb) - bwt.prefix_rank(lbRev, symb);
            size_t newLbRev = bwt.rank(lbRev, symb);
            size_t newLen   = bwt.rank(lbRev+len, symb) - newLbRev;
            auto newCursor = BiFMIndexCursor{*index, newLb, newLbRev + index->C[symb], newLen};
            return newCursor;
        }
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
        auto const& bwt = index->bwt;
        auto [rs1, prs1] = bwt.all_ranks_and_prefix_ranks(lb);
        auto [rs2, prs2] = bwt.all_ranks_and_prefix_ranks(lb+len);

        for (size_t i{0}; i < rs1.size(); ++i) {
            rs1[i] += index->C[i];
            rs2[i] += index->C[i];
        }

        auto cursors = std::array<LeftBiFMIndexCursor, Sigma>{};
        cursors[0] = LeftBiFMIndexCursor{*index, rs1[0], rs2[0] - rs1[0]};
        for (size_t i{1}; i < Sigma; ++i) {
            cursors[i] = LeftBiFMIndexCursor{*index, rs1[i], rs2[i] - rs1[i]};
        }
        return cursors;
    }

    auto extendLeft(size_t symb) const -> LeftBiFMIndexCursor {
        auto& bwt = index->bwt;

        size_t newLb    = bwt.rank(lb, symb);
        size_t newLen   = bwt.rank(lb+len, symb) - newLb;
        auto newCursor = LeftBiFMIndexCursor{*index, newLb + index->C[symb], newLen};
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
