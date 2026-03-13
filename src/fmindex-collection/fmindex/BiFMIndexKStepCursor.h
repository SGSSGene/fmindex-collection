// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "BiFMIndexKStep.h"
#include "BiFMIndexCursor.h"

namespace fmc {

template <typename Index>
struct LeftBiFMIndexKStepCursor;

template <typename Index>
struct BiFMIndexKStepCursor {
    static constexpr size_t Sigma      = Index::Sigma;
    static constexpr size_t SigmaKStep = Index::SigmaKStep;
    static constexpr size_t KStep      = Index::KStep;
    static constexpr bool   Reversed   = false;
    static constexpr size_t bitct      = std::bit_width(Sigma-1);

    constexpr bool static HasDualRank = requires(Index::String str, size_t idx) {
        { str.template all_ranks_dual_limit<bitct>(idx, idx, [](size_t, size_t, size_t, size_t, size_t) {}) };
    };

    Index const* index{};
    size_t lb;
    size_t lbRev;
    size_t len{};
    size_t steps{}; // number of extension steps taken
    BiFMIndexKStepCursor() noexcept = default;
    BiFMIndexKStepCursor(Index const& index) noexcept
        : BiFMIndexKStepCursor{index, 0, 0, index.size(), 0}
    {}
    BiFMIndexKStepCursor(Index const& index, size_t lb, size_t lbRev, size_t len, size_t steps) noexcept
        : index{&index}
        , lb{lb}
        , lbRev{lbRev}
        , len{len}
        , steps{steps}
    {}

    bool isEqual(BiFMIndexKStepCursor const& o) const noexcept {
        return lb    == o.lb
            && lbRev == o.lbRev
            && len   == o.len
            && steps == o.steps;
    }

    bool operator==(BiFMIndexKStepCursor const& _other) const noexcept {
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

    auto fetchRightBwtKStep() const -> auto const& {
        if constexpr (std::same_as<typename Index::RevBwtKStepType, std::nullptr_t>) {
            return index->bwt_kstep;
        } else {
            return index->bwtRev_kstep;
        }
    }

    auto extendLeft() const -> std::array<BiFMIndexKStepCursor, Sigma> {
        auto cursors = std::array<BiFMIndexKStepCursor, Sigma>{};
        if constexpr (HasDualRank) {
            auto& bwt = index->bwt_kstep;
            bwt.template all_ranks_dual_limit<bitct>(lb, lb+len, [&](size_t symb, size_t rs1, size_t rs2, size_t prs1, size_t prs2) {
                auto newLb    = index->C[symb] + rs1;
                auto newLen   = rs2 - rs1;
                auto newLbRev = lbRev + prs2 - prs1;
                cursors[symb] = BiFMIndexKStepCursor{*index, newLb, newLbRev, newLen, steps+1};
            });
        } else {
            for (size_t i{0}; i < Sigma; ++i) {
                cursors[i] = extendLeft(i);
            }
        }
        return cursors;
    }
    auto extendRight() const -> std::array<BiFMIndexKStepCursor, Sigma> {
        auto cursors = std::array<BiFMIndexKStepCursor, Sigma>{};
        if constexpr (HasDualRank) {
            auto& bwt = fetchRightBwtKStep();
            bwt.template all_ranks_dual_limit<bitct>(lbRev, lbRev+len, [&](size_t symb, size_t rs1, size_t rs2, size_t prs1, size_t prs2) {
                auto newLbRev = index->C[symb] + rs1;
                auto newLen   = rs2 - rs1;
                auto newLb    = lb + prs2 - prs1;
                cursors[symb] = BiFMIndexKStepCursor{*index, newLb, newLbRev, newLen, steps+1};
            });
        } else {
            for (size_t i{0}; i < Sigma; ++i) {
                cursors[i] = extendLeft(i);
            }
        }
        return cursors;
    }

    auto extendLeft(size_t symb) const -> BiFMIndexKStepCursor {
        auto& bwt = index->bwt_kstep;
        auto [lb_pr, lb_r] = bwt.template prefix_rank_and_rank_limit<bitct>(lb, symb);
        auto [rb_pr, rb_r] = bwt.template prefix_rank_and_rank_limit<bitct>(lb+len, symb);

        size_t newLb    = lb_r;
        size_t newLbRev = lbRev + rb_pr - lb_pr;
        size_t newLen   = rb_r - newLb;
        auto newCursor = BiFMIndexKStepCursor{*index, newLb + index->C[symb], newLbRev, newLen, steps+1};
        return newCursor;
    }
    auto extendRight(size_t symb) const -> BiFMIndexKStepCursor {
        auto const& bwt = fetchRightBwtKStep();

        auto [lb_pr, lb_r] = bwt.template prefix_rank_and_rank_limit<bitct>(lbRev, symb);
        auto [rb_pr, rb_r] = bwt.template prefix_rank_and_rank_limit<bitct>(lbRev+len, symb);

        size_t newLbRev = lb_r;
        size_t newLb    = lb + rb_pr - lb_pr;
        size_t newLen   = rb_r - newLbRev;
        auto newCursor = BiFMIndexKStepCursor{*index, newLb, newLbRev + index->C[symb], newLen, steps+1};
        return newCursor;
    }

    static auto convertToKSymb(std::span<size_t const, KStep> symbs) -> size_t {
        size_t kSymb{};
        for (size_t i{0}; i < KStep; ++i) {
            kSymb = kSymb*Sigma + symbs[i];
        }
        return kSymb;
    }

    /*
     * equal to extendLeft(symbs[0]).extendLeft(symbs[1])...
     */
    auto extendLeftKStep(std::span<size_t const, KStep> symbs) const -> BiFMIndexKStepCursor {
        return extendLeftKStep(convertToKSymb(symbs));
    }

    /*
     * equal to extendRight(symbs[0]).extendRight(symbs[1])...
     */
    auto extendRightKStep(std::span<size_t const, KStep> symbs) const -> BiFMIndexKStepCursor {
        return extendRightKStep(convertToKSymb(symbs));
    }

    /*
     * equal to extendLeft(symbs[0]).extendLeft(symbs[1])...
     */
    auto extendLeftKStep(size_t kSymb) const -> BiFMIndexKStepCursor {
        auto& bwt = index->bwt_kstep;
        auto [lb_pr, lb_r] = bwt.template prefix_rank_and_rank_limit<bitct*KStep>(lb, kSymb);
        auto [rb_pr, rb_r] = bwt.template prefix_rank_and_rank_limit<bitct*KStep>(lb+len, kSymb);

        size_t newLb    = lb_r;
        size_t newLbRev = lbRev + rb_pr - lb_pr;
        size_t newLen   = rb_r - newLb;
        auto newCursor = BiFMIndexKStepCursor{*index, newLb + index->C_kstep[kSymb], newLbRev, newLen, steps+KStep};
        return newCursor;
    }

    /*
     * equal to extendRight(symbs[0]).extendRight(symbs[1])...
     */
    auto extendRightKStep(size_t kSymb) const -> BiFMIndexKStepCursor {
        auto const& bwt = fetchRightBwtKStep();

        auto [lb_pr, lb_r] = bwt.template prefix_rank_and_rank_limit<bitct*KStep>(lbRev, kSymb);
        auto [rb_pr, rb_r] = bwt.template prefix_rank_and_rank_limit<bitct*KStep>(lbRev+len, kSymb);

        size_t newLbRev = lb_r;
        size_t newLb    = lb + rb_pr - lb_pr;
        size_t newLen   = rb_r - newLbRev;
        auto newCursor = BiFMIndexKStepCursor{*index, newLb, newLbRev + index->CRev_kstep[kSymb], newLen, steps+KStep};
        return newCursor;
    }


    auto extendLeftKStep() const -> std::array<BiFMIndexKStepCursor, SigmaKStep> {
        auto cursors = std::array<BiFMIndexKStepCursor, SigmaKStep>{};
        auto& bwt = index->bwt_kstep;
        if constexpr (HasDualRank) {
            bwt.template all_ranks_dual_limit<bitct*KStep>(lb, lb+len, [&](size_t symb, size_t rs1, size_t rs2, size_t prs1, size_t prs2) {
                auto newLb    = index->C_kstep[symb] + rs1;
                auto newLen   = rs2 - rs1;
                auto newLbRev = lbRev + prs2 - prs1;
                cursors[symb] = BiFMIndexKStepCursor{*index, newLb, newLbRev, newLen, steps+KStep};
            });
        } else {
            auto [rs1, prs1] = bwt.all_ranks_and_prefix_ranks(lb);
            auto [rs2, prs2] = bwt.all_ranks_and_prefix_ranks(lb+len);

            for (size_t i{0}; i < cursors.size(); ++i) {
                cursors[i] = BiFMIndexKStepCursor{*index, rs1[i] + index->C_kstep[i], lbRev + prs2[i] - prs1[i], rs2[i] - rs1[i], steps+KStep};
            }
        }
        return cursors;
    }
    auto extendRightKStep() const -> std::array<BiFMIndexKStepCursor, SigmaKStep> {
        auto cursors = std::array<BiFMIndexKStepCursor, SigmaKStep>{};
        auto& bwt = fetchRightBwtKStep();
        if constexpr (HasDualRank) {
            bwt.template all_ranks_dual_limit<bitct*KStep>(lbRev, lbRev+len, [&](size_t symb, size_t rs1, size_t rs2, size_t prs1, size_t prs2) {
                auto newLbRev = index->CRev_kstep[symb] + rs1;
                auto newLen   = rs2 - rs1;
                auto newLb    = lb + prs2 - prs1;
                cursors[symb] = BiFMIndexKStepCursor{*index, newLb, newLbRev, newLen, steps+KStep};
            });
        } else {
            auto [rs1, prs1] = bwt.all_ranks_and_prefix_ranks(lbRev);
            auto [rs2, prs2] = bwt.all_ranks_and_prefix_ranks(lbRev+len);

            for (size_t i{0}; i < cursors.size(); ++i) {
                cursors[i] = BiFMIndexKStepCursor{*index, lb + prs2[i] - prs1[i], rs1[i] + index->CRev_kstep[i], rs2[i] - rs1[i], steps+KStep};
            }
        }
        return cursors;
    }

    // This requires that all rows have the same BWT entry (or only a single one is available)
    // - must have at least marked a single row
    // - all rows must have the same 'bwt' symbol
    auto extendLeftBySymbolKStep(size_t kSymb) const -> BiFMIndexKStepCursor {
        auto& bwt = index->bwt_kstep;

        assert(count() > 0);
        assert([&]() {
            for (size_t i{lb}; i < lb + len; ++i) {
                if (kSymb != bwt.symbol(i)) {
                    return false;
                }
            }
            return true;
        }());

        size_t newLb    = bwt.rank(lb, kSymb);
        size_t newLbRev = lbRev;
        size_t newLen   = len;
        auto newCursor  = BiFMIndexKStepCursor{*index, newLb + index->C_kstep[kSymb], newLbRev, newLen, steps+1};
        return newCursor;
    }

    // see extendLeftBySymbol
    auto extendRightBySymbolKStep(size_t kSymb) const -> BiFMIndexKStepCursor {
        auto const& bwt = fetchRightBwtKStep();
        size_t newLb    = lb;
        size_t newLbRev = bwt.rank(lbRev, kSymb);
        size_t newLen   = len;
        auto newCursor = BiFMIndexKStepCursor{*index, newLb, newLbRev + index->CRev_kstep[kSymb], newLen, steps+1};
        return newCursor;
    }


    // This requires that all rows have the same BWT entry (or only a single one is available)
    // - must have at least marked a single row
    // - all rows must have the same 'bwt' symbol
    auto extendLeftBySymbolKStep() const -> std::tuple<size_t, BiFMIndexKStepCursor> {
        auto& bwt = index->bwt_kstep;
        auto symb = bwt.symbol(lb);
        return {symb, extendLeftBySymbolKStep(symb)};
    }

    // see extendLeftBySymbol
    auto extendRightBySymbolKStep() const -> std::tuple<size_t, BiFMIndexKStepCursor> {
        auto const& bwt = fetchRightBwtKStep();
        auto symb = bwt.symbol(lbRev);
        return {symb, extendRightBySymbolKStep(symb)};
    }



    auto symbolLeft() const -> size_t {
        return index->bwt_kstep.template symbol_limit<bitct>(lb);
    }
    auto symbolRight() const -> size_t {
        auto const& bwt = fetchRightBwtKStep();
        return bwt.template symbol_limit<bitct>(lbRev);
    }

    auto symbolLeft(size_t idx) const -> size_t {
        return index->bwt_kstep.template symbol_limit<bitct>(lb+idx);
    }
    auto symbolRight(size_t idx) const -> size_t {
        auto const& bwt = fetchRightBwtKStep();
        return bwt.template symbol_limit<bitct>(lbRev+idx);
    }


    auto symbolLeftKStep() const -> size_t {
        return index->bwt_kstep.symbol(lb);
    }
    auto symbolRightKStep() const -> size_t {
        auto const& kbwt = fetchRightBwtKStep();
        return kbwt.symbol(lbRev);
    }

};

template <typename Index>
auto begin(BiFMIndexKStepCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb};
}
template <typename Index>
auto end(BiFMIndexKStepCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb + _cursor.len};
}

template <typename Index>
struct LeftBiFMIndexKStepCursor {
    static constexpr size_t Sigma      = Index::Sigma;
    static constexpr size_t SigmaKStep = Index::SigmaKStep;
    static constexpr size_t KStep      = Index::KStep;
    static constexpr bool   Reversed = false;

    Index const* index;
    size_t lb;
    size_t len;
    size_t steps;
    LeftBiFMIndexKStepCursor(BiFMIndexKStepCursor<Index> const& _other)
        : index{_other.index}
        , lb{_other.lb}
        , len{_other.len}
        , steps{_other.steps}
    {}
    LeftBiFMIndexKStepCursor()
        : index{nullptr}
    {}
    LeftBiFMIndexKStepCursor(Index const& index)
        : LeftBiFMIndexKStepCursor{index, 0, index.size(), 0}
    {}
    LeftBiFMIndexKStepCursor(Index const& index, size_t lb, size_t len, size_t steps)
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

    auto extendLeft() const -> std::array<LeftBiFMIndexKStepCursor, Sigma> {
        auto cursors = std::array<LeftBiFMIndexKStepCursor, Sigma>{};

        for (size_t i{0}; i < Sigma; ++i) {
            cursors[i] = extendLeft(i);
        }
        return cursors;
    }

    auto extendLeft(size_t symb) const -> LeftBiFMIndexKStepCursor {
        auto& bwt = index->bwt_kstep;

        size_t newLb    = bwt.rank(lb, symb);
        size_t newLen   = bwt.rank(lb+len, symb) - newLb;
        auto newCursor = LeftBiFMIndexKStepCursor{*index, newLb + index->C[symb], newLen, steps+1};
        return newCursor;
    }

    auto extendLeftKStep(std::span<size_t const, KStep> symbs) const -> LeftBiFMIndexKStepCursor {
        size_t symb_pr{};
        for (size_t i{0}; i < KStep; ++i) {
            symb_pr = symb_pr*Sigma + symbs[i];
        }

        auto& bwt = index->bwt_kstep;
        size_t newLb    = bwt.rank(lb, symb_pr);
        size_t newLen   = bwt.rank(lb+len, symb_pr) - newLb;
        auto newCursor = LeftBiFMIndexKStepCursor{*index, newLb + index->C_kstep[symb_pr], newLen, steps+KStep};
        return newCursor;
    }
};

template <typename Index>
auto begin(LeftBiFMIndexKStepCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb};
}
template <typename Index>
auto end(LeftBiFMIndexKStepCursor<Index> const& _cursor) {
    return IntIterator{_cursor.lb + _cursor.len};
}

}

namespace std {

template <typename index_t>
struct hash<fmc::BiFMIndexKStepCursor<index_t>> {
    auto operator()(fmc::BiFMIndexKStepCursor<index_t> const& cursor) const -> size_t {
        return hash<size_t>()(cursor.lb)
            ^ hash<size_t>()(cursor.len);
    }
};

}
