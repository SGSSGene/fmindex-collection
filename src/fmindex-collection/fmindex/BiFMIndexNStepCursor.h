// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "BiFMIndexNStep.h"
#include "BiFMIndexCursor.h"

namespace fmc {

template <typename Index>
struct LeftBiFMIndexNStepCursor;

template <typename Index>
struct BiFMIndexNStepCursor {
    static constexpr size_t Sigma      = Index::Sigma;
    static constexpr size_t SigmaNStep = Index::SigmaNStep;
    static constexpr size_t NStep      = Index::NStep;
    static constexpr bool   Reversed = false;

    Index const* index{};
    size_t lb;
    size_t lbRev;
    size_t len{};
    size_t steps{}; // number of extension steps taken
    BiFMIndexNStepCursor() noexcept = default;
    BiFMIndexNStepCursor(Index const& index) noexcept
        : BiFMIndexNStepCursor{index, 0, 0, index.size(), 0}
    {}
    BiFMIndexNStepCursor(Index const& index, size_t lb, size_t lbRev, size_t len, size_t steps) noexcept
        : index{&index}
        , lb{lb}
        , lbRev{lbRev}
        , len{len}
        , steps{steps}
    {}

    bool operator==(BiFMIndexNStepCursor const& _other) const noexcept {
        return lb == _other.lb
               && len == _other.len;
    }

    operator BiFMIndexCursor<Index>() const {
        return BiFMIndexCursor{*index, lb, lbRev, len, steps};
    }
    operator LeftBiFMIndexCursor<Index>() const {
        return LeftBiFMIndexCursor{*index, lb, len, steps};
    }

    bool empty() const {
        return len == 0;
    }
    size_t count() const {
        return len;
    }
    auto extendLeft() const -> std::array<BiFMIndexNStepCursor, Sigma> {
        auto const& bwt = index->bwt;
        auto [rs1, prs1] = bwt.all_ranks_and_prefix_ranks(lb);
        auto [rs2, prs2] = bwt.all_ranks_and_prefix_ranks(lb+len);

        for (size_t i{0}; i < rs1.size(); ++i) {
            rs1[i] += index->C[i];
            rs2[i] += index->C[i];
        }

        auto cursors = std::array<BiFMIndexNStepCursor, Sigma>{};
        for (size_t i{0}; i < cursors.size(); ++i) {
            cursors[i] = BiFMIndexNStepCursor{*index, rs1[i], lbRev + prs2[i] - prs1[i], rs2[i] - rs1[i], steps+1};
        }
        return cursors;
    }

    auto extendRight() const -> std::array<BiFMIndexNStepCursor, Sigma> {
        // Reuse bwt if bwtrev is not available
        if constexpr (std::same_as<typename Index::RevBwtType, std::nullptr_t>) {
            auto const& bwt = index->bwt;
            auto [rs1, prs1] = bwt.all_ranks_and_prefix_ranks(lbRev);
            auto [rs2, prs2] = bwt.all_ranks_and_prefix_ranks(lbRev+len);

            for (size_t i{0}; i < rs1.size(); ++i) {
                rs1[i] += index->C[i];
                rs2[i] += index->C[i];
            }

            auto cursors = std::array<BiFMIndexNStepCursor, Sigma>{};
            for (size_t i{0}; i < cursors.size(); ++i) {
                cursors[i] = BiFMIndexNStepCursor{*index, lb + prs2[i] - prs1[i], rs1[i], rs2[i] - rs1[i], steps+1};
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

            auto cursors = std::array<BiFMIndexNStepCursor, Sigma>{};
            for (size_t i{0}; i < cursors.size(); ++i) {
                cursors[i] = BiFMIndexNStepCursor{*index, lb + prs2[i] - prs1[i], rs1[i], rs2[i] - rs1[i], steps+1};
            }
            return cursors;
        }
    }

    void prefetchLeft() const {
    }
    void prefetchRight() const {
    }

    auto extendLeft(size_t symb) const -> BiFMIndexNStepCursor {
        auto& bwt = index->bwt;
        size_t newLb    = bwt.rank(lb, symb);
        size_t newLbRev = lbRev + bwt.prefix_rank(lb+len, symb) - bwt.prefix_rank(lb, symb);
        size_t newLen   = bwt.rank(lb+len, symb) - newLb;
        auto newCursor = BiFMIndexNStepCursor{*index, newLb + index->C[symb], newLbRev, newLen, steps+1};
        return newCursor;
    }
    auto extendRight(size_t symb) const -> BiFMIndexNStepCursor {
        if constexpr (std::same_as<typename Index::RevBwtType, std::nullptr_t>) {
            auto const& bwt = index->bwt;
            size_t newLb    = lb + bwt.prefix_rank(lbRev+len, symb) - bwt.prefix_rank(lbRev, symb);
            size_t newLbRev = bwt.rank(lbRev, symb);
            size_t newLen   = bwt.rank(lbRev+len, symb) - newLbRev;
            auto newCursor = BiFMIndexNStepCursor{*index, newLb, newLbRev + index->C[symb], newLen, steps+1};
            return newCursor;
        } else {
            auto const& bwt = index->bwtRev;
            size_t newLb    = lb + bwt.prefix_rank(lbRev+len, symb) - bwt.prefix_rank(lbRev, symb);
            size_t newLbRev = bwt.rank(lbRev, symb);
            size_t newLen   = bwt.rank(lbRev+len, symb) - newLbRev;
            auto newCursor = BiFMIndexNStepCursor{*index, newLb, newLbRev + index->C[symb], newLen, steps+1};
            return newCursor;
        }
    }

    /*
     * equal to extendLeft(symbs[0]).extendLeft(symbs[1])...
     */
    auto extendLeftNStep(std::span<size_t, NStep> symbs) const -> BiFMIndexNStepCursor {
        size_t symb_pr{};
        for (size_t i{0}; i < NStep; ++i) {
            symb_pr = symb_pr*Sigma + symbs[i];
        }

        auto& bwt = index->bwt_nstep;
        size_t newLb    = bwt.rank(lb, symb_pr);
        size_t newLbRev = lbRev + bwt.prefix_rank(lb+len, symb_pr) - bwt.prefix_rank(lb, symb_pr);
        size_t newLen   = bwt.rank(lb+len, symb_pr) - newLb;
        auto newCursor = BiFMIndexNStepCursor{*index, newLb + index->C_nstep[symb_pr], newLbRev, newLen, steps+NStep};
        return newCursor;
    }

    /*
     * equal to extendRight(symbs[0]).extendRight(symbs[1])...
     */
    auto extendRightNStep(std::span<size_t, NStep> symbs) const -> BiFMIndexNStepCursor {
        size_t symb_pr{};
        for (size_t i{0}; i < NStep; ++i) {
            symb_pr = symb_pr*Sigma + symbs[i];
        }
        if constexpr (std::same_as<typename Index::RevBwtNStepType, std::nullptr_t>) {
            auto const& bwt = index->bwt_nstep;
            size_t newLb    = lb + bwt.prefix_rank(lbRev+len, symb_pr) - bwt.prefix_rank(lbRev, symb_pr);
            size_t newLbRev = bwt.rank(lbRev, symb_pr);
            size_t newLen   = bwt.rank(lbRev+len, symb_pr) - newLbRev;
            auto newCursor = BiFMIndexNStepCursor{*index, newLb, newLbRev + index->C_nstep[symb_pr], newLen, steps+NStep};
            return newCursor;
        } else {
            auto const& bwt = index->bwtRev_nstep;
            size_t newLb    = lb + bwt.prefix_rank(lbRev+len, symb_pr) - bwt.prefix_rank(lbRev, symb_pr);
            size_t newLbRev = bwt.rank(lbRev, symb_pr);
            size_t newLen   = bwt.rank(lbRev+len, symb_pr) - newLbRev;
            auto newCursor = BiFMIndexNStepCursor{*index, newLb, newLbRev + index->CRev_nstep[symb_pr], newLen, steps+NStep};
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

    auto symbolLeftNStep() const -> size_t {
        return index->bwt_nstep.symbol(lb);
    }
    auto symbolRightNStep() const -> size_t {
        // Reuse bwt if bwtrev is not available
        if constexpr (std::same_as<typename Index::RevBwtNStepType, std::nullptr_t>) {
            return index->bwt_nstep.symbol(lbRev);
        } else {
            return index->bwtRev_nstep.symbol(lbRev);
        }
    }

};

template <typename Index>
auto begin(BiFMIndexNStepCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb};
}
template <typename Index>
auto end(BiFMIndexNStepCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb + _cursor.len};
}

template <typename Index>
struct LeftBiFMIndexNStepCursor {
    static constexpr size_t Sigma    = Index::Sigma;
    static constexpr bool   Reversed = false;

    Index const* index;
    size_t lb;
    size_t len;
    size_t steps;
    LeftBiFMIndexNStepCursor(BiFMIndexNStepCursor<Index> const& _other)
        : index{_other.index}
        , lb{_other.lb}
        , len{_other.len}
        , steps{_other.steps}
    {}
    LeftBiFMIndexNStepCursor()
        : index{nullptr}
    {}
    LeftBiFMIndexNStepCursor(Index const& index)
        : LeftBiFMIndexNStepCursor{index, 0, index.size(), 0}
    {}
    LeftBiFMIndexNStepCursor(Index const& index, size_t lb, size_t len, size_t steps)
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

    operator LeftBiFMIndexCursor<Index>() const {
        return LeftBiFMIndexCursor{*index, lb, len, steps};
    }

    auto extendLeft() const -> std::array<LeftBiFMIndexNStepCursor, Sigma> {
        auto const& bwt = index->bwt;
        auto [rs1, prs1] = bwt.all_ranks_and_prefix_ranks(lb);
        auto [rs2, prs2] = bwt.all_ranks_and_prefix_ranks(lb+len);

        for (size_t i{0}; i < rs1.size(); ++i) {
            rs1[i] += index->C[i];
            rs2[i] += index->C[i];
        }

        auto cursors = std::array<LeftBiFMIndexNStepCursor, Sigma>{};
        cursors[0] = LeftBiFMIndexNStepCursor{*index, rs1[0], rs2[0] - rs1[0], steps+1};
        for (size_t i{1}; i < Sigma; ++i) {
            cursors[i] = LeftBiFMIndexNStepCursor{*index, rs1[i], rs2[i] - rs1[i], steps+1};
        }
        return cursors;
    }

    auto extendLeft(size_t symb) const -> LeftBiFMIndexNStepCursor {
        auto& bwt = index->bwt;

        size_t newLb    = bwt.rank(lb, symb);
        size_t newLen   = bwt.rank(lb+len, symb) - newLb;
        auto newCursor = LeftBiFMIndexNStepCursor{*index, newLb + index->C[symb], newLen, steps+1};
        return newCursor;
    }
};

template <typename Index>
auto begin(LeftBiFMIndexNStepCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb};
}
template <typename Index>
auto end(LeftBiFMIndexNStepCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb + _cursor.len};
}

}

namespace std {

template <typename index_t>
struct hash<fmc::BiFMIndexNStepCursor<index_t>> {
    auto operator()(fmc::BiFMIndexNStepCursor<index_t> const& cursor) const -> size_t {
        return hash<size_t>()(cursor.lb)
            ^ hash<size_t>()(cursor.len);
    }
};

}
