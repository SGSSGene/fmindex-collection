// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "BinaryMirroredBiFMIndex.h"

namespace fmc {

template <typename Index>
struct LeftBinaryMirroredBiFMIndexCursor;

template <typename Index>
struct BinaryMirroredBiFMIndexCursor {
    static constexpr size_t Sigma    = Index::Sigma;
    static constexpr bool   Reversed = false;

    Index const* index{};
    size_t lb;
    size_t lbRev;
    size_t len{};
    BinaryMirroredBiFMIndexCursor() noexcept = default;
    BinaryMirroredBiFMIndexCursor(Index const& index) noexcept
        : BinaryMirroredBiFMIndexCursor{index, 0, 0, index.size()}
    {}
    BinaryMirroredBiFMIndexCursor(Index const& index, size_t lb, size_t lbRev, size_t len) noexcept
        : index{&index}
        , lb{lb}
        , lbRev{lbRev}
        , len{len}
    {}

    bool operator==(BinaryMirroredBiFMIndexCursor const& _other) const noexcept {
        return lb == _other.lb
               && len == _other.len;
    }
    bool empty() const {
        return len == 0;
    }
    size_t count() const {
        return len;
    }

private:
    // Helper function, since Bitvectors don't have prefix rank operations and not symbol specific rank
    auto ranks_and_prefixes(size_t i) const -> std::tuple<std::array<size_t, 2>, std::array<size_t, 2>> {
        auto r = index->bwt.rank(i);

        auto rs  = std::array<size_t, 2>{i-r, r};
        auto prs = std::array<size_t, 2>{rs[0], i};
        rs[0] += index->C[0];
        rs[1] += index->C[1];
        return {rs, prs};
    };

    auto rank(size_t i, size_t symb) const -> size_t {
        size_t r = index->bwt.rank(i, symb);
        return (symb*r) + (1-symb)*(i-r);
    }

public:
    auto extendLeft() const -> std::array<BinaryMirroredBiFMIndexCursor, Sigma> {
        auto [rs1, prs1] = ranks_and_prefixes(lb);
        auto [rs2, prs2] = ranks_and_prefixes(lb+len);

        auto cursors = std::array<BinaryMirroredBiFMIndexCursor, Sigma>{};
        cursors[0] = BinaryMirroredBiFMIndexCursor{*index, rs1[0], lbRev, rs2[0] - rs1[0]};
        for (size_t i{1}; i < Sigma; ++i) {
            cursors[i] = BinaryMirroredBiFMIndexCursor{*index, rs1[i], lbRev + prs2[i-1] - prs1[i-1], rs2[i] - rs1[i]};
        }
        return cursors;
    }

    auto extendRight() const -> std::array<BinaryMirroredBiFMIndexCursor, Sigma> {
        auto [rs1, prs1] = ranks_and_prefixes(lbRev);
        auto [rs2, prs2] = ranks_and_prefixes(lbRev+len);

        auto cursors = std::array<BinaryMirroredBiFMIndexCursor, Sigma>{};
        cursors[0] = BinaryMirroredBiFMIndexCursor{*index, lb, rs1[0], rs2[0] - rs1[0]};
        cursors[0].prefetchRight();
        for (size_t i{1}; i < Sigma; ++i) {
            cursors[i] = BinaryMirroredBiFMIndexCursor{*index, lb + prs2[i-1] - prs1[i-1], rs1[i], rs2[i] - rs1[i]};
        }
        return cursors;
    }
    void prefetchLeft() const {
    }
    void prefetchRight() const {
    }

    auto extendLeft(size_t symb) const -> BinaryMirroredBiFMIndexCursor {
        size_t newLb    = rank(lb, symb);
        size_t newLbRev = lbRev + [&]() -> size_t {
            if (symb == 0) return {};
            return index->bwt.rank(lb+len) - index->bwt.rank(lb);
        }();
        size_t newLen   = rank(lb+len, symb) - newLb;
        auto newCursor = BinaryMirroredBiFMIndexCursor{*index, newLb + index->C[symb], newLbRev, newLen};
        return newCursor;
    }
    auto extendRight(size_t symb) const -> BinaryMirroredBiFMIndexCursor {
        size_t newLb    = lb + [&]() -> size_t {
            if (symb == 0) return {};
            return index->bwt.rank(lbRev+len) - index->bwt.rank(lbRev);
        }();
        size_t newLbRev = rank(lbRev, symb);
        size_t newLen   = rank(lbRev+len, symb) - newLbRev;
        auto newCursor = BinaryMirroredBiFMIndexCursor{*index, newLb, newLbRev + index->C[symb], newLen};
        return newCursor;
    }
};

template <typename Index>
auto begin(BinaryMirroredBiFMIndexCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb};
}
template <typename Index>
auto end(BinaryMirroredBiFMIndexCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb + _cursor.len};
}

template <typename Index>
struct LeftBinaryMirroredBiFMIndexCursor {
    static constexpr size_t Sigma    = Index::Sigma;
    static constexpr bool   Reversed = false;

    Index const* index;
    size_t lb;
    size_t len;
    LeftBinaryMirroredBiFMIndexCursor(BinaryMirroredBiFMIndexCursor<Index> const& _other)
        : index{_other.index}
        , lb{_other.lb}
        , len{_other.len}
    {}
    LeftBinaryMirroredBiFMIndexCursor()
        : index{nullptr}
    {}
    LeftBinaryMirroredBiFMIndexCursor(Index const& index)
        : LeftBinaryMirroredBiFMIndexCursor{index, 0, index.size()}
    {}
    LeftBinaryMirroredBiFMIndexCursor(Index const& index, size_t lb, size_t len)
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

private:
    // Helper function, since Bitvectors don't have prefix rank operations and not symbol specific rank
    auto ranks_and_prefixes(size_t i) const -> std::tuple<std::array<size_t, 2>, std::array<size_t, 2>> {
        auto r = index->bwt.rank(i);

        auto rs  = std::array<size_t, 2>{i-r, r};
        auto prs = std::array<size_t, 2>{rs[0], i};
        rs[0] += index->C[0];
        rs[1] += index->C[1];
        return {rs, prs};
    };

    auto rank(size_t i, size_t symb) const -> size_t {
        size_t r = index->bwt.rank(i, symb);
        return (symb*r) + (1-symb)*(i-r);
    }

public:

    auto extendLeft() const -> std::array<LeftBinaryMirroredBiFMIndexCursor, Sigma> {
        auto [rs1, prs1] = ranks_and_prefixes(lb);
        auto [rs2, prs2] = ranks_and_prefixes(lb+len);

        auto cursors = std::array<LeftBinaryMirroredBiFMIndexCursor, Sigma>{};
        cursors[0] = LeftBinaryMirroredBiFMIndexCursor{*index, rs1[0], rs2[0] - rs1[0]};
        for (size_t i{1}; i < Sigma; ++i) {
            cursors[i] = LeftBinaryMirroredBiFMIndexCursor{*index, rs1[i], rs2[i] - rs1[i]};
        }
        return cursors;
    }

    auto extendLeft(size_t symb) const -> LeftBinaryMirroredBiFMIndexCursor {
        size_t newLb    = rank(lb, symb);
        size_t newLen   = rank(lb+len, symb) - newLb;

        auto newCursor = LeftBinaryMirroredBiFMIndexCursor{*index, newLb + index->C[symb], newLen};
        return newCursor;
    }
};

template <typename Index>
auto begin(LeftBinaryMirroredBiFMIndexCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb};
}
template <typename Index>
auto end(LeftBinaryMirroredBiFMIndexCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb + _cursor.len};
}

}

namespace std {

template <typename index_t>
struct hash<fmc::BinaryMirroredBiFMIndexCursor<index_t>> {
    auto operator()(fmc::BinaryMirroredBiFMIndexCursor<index_t> const& cursor) const -> size_t {
        return hash<size_t>()(cursor.lb)
            ^ hash<size_t>()(cursor.len);
    }
};

}
