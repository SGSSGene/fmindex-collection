// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../fmindex/BiFMIndex.h"
#include "../fmindex/BiFMIndexCursor.h"
#include "../fmindex/LinearFMIndex.h"
#include "../fmindex/LinearFMIndexCursor.h"
#include "../search_scheme/Scheme.h"
#include "Restore.h"
#include "SelectCursor.h"
#include "concepts.h"


// Double Index using a bidirectional fm-index and a linear fm index
namespace fmc::search_double_index2 {

template <bool Edit, typename Index, typename QIndex, typename search_t, typename delegate_t>
struct Search {

    using cursor_t  = BiFMIndexCursor<Index>;
    using qcursor_t = LinearFMIndexCursor<QIndex>;

    Index  const& index;
    QIndex const& queries;

    decltype(search_t::pi) const& pi;
    decltype(search_t::l) const& l;
    decltype(search_t::u) const& u;

    size_t threshold;
    size_t optMode;

    delegate_t const& delegate;

    struct Side {
        uint8_t lastRank{};
        uint8_t lastQRank{};
    };
    mutable std::array<Side, 2> side;
    mutable size_t e{};
    mutable size_t part{};

    Search(Index const& _index, QIndex const& _queries, search_t const& _search, size_t _threshold, size_t _optMode, delegate_t const& _delegate) noexcept
        : index {_index}
        , queries{_queries}
        , pi{_search.pi}
        , l{_search.l}
        , u{_search.u}
        , threshold{_threshold}
        , optMode{_optMode}
        , delegate  {_delegate}
    {
        auto cur       = cursor_t{index};
        auto qcur      = qcursor_t{queries};
        searchPart<'M', 'M'>(cur, qcur);
    }

    template <char LInfo, char RInfo>
    void searchPart(cursor_t const& cur, qcursor_t const& qcur) const noexcept {
        if (cur.count() == 0 || qcur.count() == 0) {
            return;
        }

        if (part == pi.size()) {
            if (l[part-1] <= e and e <= u[part-1]) {
                if constexpr (!Edit || ((LInfo == 'M' or LInfo == 'I') and (RInfo == 'M' or RInfo == 'I'))) {
                    delegate(cur, qcur, e);
                }
            }
            return;
        }
        if (e > u[part]) {
            return;
        }

        if (optMode & 0x04) {
            if (cur.count() == 1 && qcur.count() == 1) {
                if (part == 0 || pi[part-1] < pi[part]) {
                    searchPartDirSingleQuerySingleIndex<LInfo, RInfo, true>(cur, qcur);
                } else {
                    searchPartDirSingleQuerySingleIndex<LInfo, RInfo, false>(cur, qcur);
                }
                return;

            }
        }
        if (optMode & 0x01) {
            if (qcur.count() <= threshold) {
                if (qcur.count() == 1) {
                    if (part == 0 || pi[part-1] < pi[part]) {
                        searchPartDirSingleQuery<LInfo, RInfo, true>(cur, qcur);
                    } else {
                        searchPartDirSingleQuery<LInfo, RInfo, false>(cur, qcur);
                    }
                    return;
                }

                auto short_qcur = qcur;
                short_qcur.len = 1;
                for (size_t i{}; i < qcur.count(); ++i) {
                    short_qcur.lb = qcur.lb + i;
                    searchPart<LInfo, RInfo>(cur, short_qcur);
                }
                return;
            }
        }
        if (optMode & 0x02) {
            if (cur.count() == 1) {
                if (part == 0 || pi[part-1] < pi[part]) {
                    searchPartDirSingleIndex<LInfo, RInfo, true>(cur, qcur);
                } else {
                    searchPartDirSingleIndex<LInfo, RInfo, false>(cur, qcur);
                }
                return;
            }
        }

        if (part == 0 || pi[part-1] < pi[part]) {
            searchPartDir<LInfo, RInfo, true>(cur, qcur);
        } else {
            searchPartDir<LInfo, RInfo, false>(cur, qcur);
        }
    }

    template <char LInfo, char RInfo, bool Right>
    void searchPartDir(cursor_t const& cur, qcursor_t const& qcur) const noexcept {
        static constexpr char TInfo = Right ? RInfo : LInfo;

        constexpr bool Deletion     = (TInfo != 'S' && TInfo != 'I') && Edit;
        constexpr bool Insertion    = (TInfo != 'S' && TInfo != 'D') && Edit;

        constexpr char OnMatchL      = Right ? LInfo : 'M';
        constexpr char OnMatchR      = Right ? 'M'   : RInfo;
        constexpr char OnSubstituteL = Right ? LInfo : 'S';
        constexpr char OnSubstituteR = Right ? 'S'   : RInfo;
        constexpr char OnDeletionL   = Right ? LInfo : 'D';
        constexpr char OnDeletionR   = Right ? 'D'   : RInfo;
        constexpr char OnInsertionL  = Right ? LInfo : 'I';
        constexpr char OnInsertionR  = Right ? 'I'   : RInfo;

        auto lastRank  = side[Right].lastRank;
        auto lastQRank = side[Right].lastQRank;

        auto matchAllowed    = [lastRank, lastQRank](auto nextSymb) {
                                return (TInfo != 'I' or nextSymb != lastQRank)
                                   and (TInfo != 'D'  or nextSymb != lastRank);
                                };
        bool mismatchAllowed = l[part] <= e+1 and e+1 <= u[part];

        auto qcursors = qcur.extendLeft();
        auto cursors = [&]() {
            if constexpr (Right) {
                return cur.extendRight();
            } else {
                return cur.extendLeft();
            }
        }();

        // search matches
        if (l[part] <= e and e <= u[part]) {
            auto sp = Restore{part, part+1};
            auto s1 = Restore{side[Right].lastRank};
            auto s2 = Restore{side[Right].lastQRank};

            for (size_t symb{1}; symb < std::min(Index::Sigma, QIndex::Sigma); ++symb) {
                if (matchAllowed(symb)) {
                    side[Right].lastRank = symb;
                    side[Right].lastQRank = symb;
                    searchPart<OnMatchL, OnMatchR>(cursors[symb], qcursors[symb]);
                }
            }
        }
        if (mismatchAllowed) {
            auto se = Restore{e, e+1};
            // search substitutes
            {
                auto s1 = Restore{side[Right].lastRank};
                auto s2 = Restore{side[Right].lastQRank};
                auto sp = Restore{part, part+1};
                for (size_t qsymb{1}; qsymb < QIndex::Sigma; ++qsymb) {
                    if (qcursors[qsymb].count() == 0) continue;
                    side[Right].lastQRank = qsymb;
                    for (size_t symb{1}; symb < Index::Sigma; ++symb) {
                        if (symb == qsymb) continue;
                        side[Right].lastRank = symb;
                        searchPart<OnSubstituteL, OnSubstituteR>(cursors[symb], qcursors[qsymb]);
                    }
                }
            }
            // search deletion
            if constexpr (Deletion) {
                auto s1 = Restore{side[Right].lastRank};
                for (size_t symb{1}; symb < Index::Sigma; ++symb) {
                    if (TInfo != 'M' || lastQRank != symb) {
                        side[Right].lastRank = symb;
                        searchPart<OnDeletionL, OnDeletionR>(cursors[symb], qcur);
                    }
                }
            }

            // search insertion
            if constexpr (Insertion) {
                auto s1 = Restore{side[Right].lastQRank};
                auto sp = Restore{part, part+1};
                for (size_t qsymb{1}; qsymb < QIndex::Sigma; ++qsymb) {
                    if (TInfo != 'M' || lastQRank != qsymb) {
                        side[Right].lastQRank = qsymb;
                        searchPart<OnInsertionL, OnInsertionR>(cur, qcursors[qsymb]);
                    }
                }
            }
        }
    }
#if 1
    template <char LInfo, char RInfo, bool Right>
    void searchPartDirSingleQuerySingleIndex(cursor_t const& cur, qcursor_t const& qcur) const noexcept {
        static constexpr char TInfo = Right ? RInfo : LInfo;

        constexpr bool Deletion     = (TInfo != 'S' && TInfo != 'I') && Edit;
        constexpr bool Insertion    = (TInfo != 'S' && TInfo != 'D') && Edit;

        constexpr char OnMatchL      = Right ? LInfo : 'M';
        constexpr char OnMatchR      = Right ? 'M'   : RInfo;
        constexpr char OnSubstituteL = Right ? LInfo : 'S';
        constexpr char OnSubstituteR = Right ? 'S'   : RInfo;
        constexpr char OnDeletionL   = Right ? LInfo : 'D';
        constexpr char OnDeletionR   = Right ? 'D'   : RInfo;
        constexpr char OnInsertionL  = Right ? LInfo : 'I';
        constexpr char OnInsertionR  = Right ? 'I'   : RInfo;

        auto lastRank  = side[Right].lastRank;
        auto lastQRank = side[Right].lastQRank;

        auto [curISymb, icursorNext] = [&]() -> std::tuple<size_t, cursor_t> {
            if constexpr (Right) {
                auto symb = cur.symbolRight();
                auto cur_  = cur.extendRight(symb);
                return {symb, cur_};
            } else {
                auto symb = cur.symbolLeft();
                auto cur_  = cur.extendLeft(symb);
                return {symb, cur_};
            }
        }();

        auto curQSymb = queries.columns[qcur.col].bwt.symbol(qcur.lb);
        auto qcursorNext = qcur.extendLeft(curQSymb);

        bool matchAllowed    = l[part] <= e and e <= u[part]
                                and (TInfo != 'I' or curQSymb != lastQRank)
                                and (TInfo != 'D'  or curQSymb != lastRank);

        bool mismatchAllowed = l[part] <= e+1 and e+1 <= u[part];


        // only insertions matter
        if constexpr (Insertion) {
            if (mismatchAllowed) {
                if (TInfo != 'M' || lastQRank != curQSymb) {
                    auto s1 = Restore{side[Right].lastQRank, curQSymb};
                    auto sp = Restore{part, part+1};
                    auto se = Restore{e, e+1};
                    searchPart<OnInsertionL, OnInsertionR>(cur, qcursorNext);
                }
            }
        }

        if (curISymb == 0) return;

        // match
        if (curQSymb == curISymb && matchAllowed) {
            auto sp = Restore{part, part+1};
            auto s1 = Restore{side[Right].lastRank, curISymb};
            auto s2 = Restore{side[Right].lastQRank, curQSymb};
            searchPart<OnMatchL, OnMatchR>(icursorNext, qcursorNext);
        }
        if (mismatchAllowed) {
            auto se = Restore{e, e+1};
            // search substitutes
            if (curQSymb != curISymb)
            {
                auto s1 = Restore{side[Right].lastRank, curISymb};
                auto s2 = Restore{side[Right].lastQRank, curQSymb};
                auto sp = Restore{part, part+1};
                searchPart<OnSubstituteL, OnSubstituteR>(icursorNext, qcursorNext);
            }
            // search deletion
            if constexpr (Deletion) {
                if (TInfo != 'M' || lastQRank != curISymb) {
                    auto s1 = Restore{side[Right].lastRank, curISymb};
                    searchPart<OnDeletionL, OnDeletionR>(icursorNext, qcur);
                }
            }
        }
    }
#endif
    #if 1
    template <char LInfo, char RInfo, bool Right>
    void searchPartDirSingleQuery(cursor_t const& cur, qcursor_t const& qcur) const noexcept {
        static constexpr char TInfo = Right ? RInfo : LInfo;

        constexpr bool Deletion     = (TInfo != 'S' && TInfo != 'I') && Edit;
        constexpr bool Insertion    = (TInfo != 'S' && TInfo != 'D') && Edit;

        constexpr char OnMatchL      = Right ? LInfo : 'M';
        constexpr char OnMatchR      = Right ? 'M'   : RInfo;
        constexpr char OnSubstituteL = Right ? LInfo : 'S';
        constexpr char OnSubstituteR = Right ? 'S'   : RInfo;
        constexpr char OnDeletionL   = Right ? LInfo : 'D';
        constexpr char OnDeletionR   = Right ? 'D'   : RInfo;
        constexpr char OnInsertionL  = Right ? LInfo : 'I';
        constexpr char OnInsertionR  = Right ? 'I'   : RInfo;


        auto curQSymb = queries.columns[qcur.col].bwt.symbol(qcur.lb);
        auto qcursorNext = qcur.extendLeft(curQSymb);


        bool matchAllowed    = l[part] <= e and e <= u[part]
                               and ((TInfo != 'I' or curQSymb != side[Right].lastQRank)
                               and (TInfo != 'D' or curQSymb != side[Right].lastRank));
        bool mismatchAllowed = l[part] <= e+1 and e+1 <= u[part];

        // search matches
        if (matchAllowed) {
            auto cursor = [&]() {
                if constexpr (Right) {
                    return cur.extendRight(curQSymb);
                } else {
                    return cur.extendLeft(curQSymb);
                }
            }();

            auto sp = Restore{part, part+1};
            auto s1 = Restore{side[Right].lastRank, curQSymb};
            auto s2 = Restore{side[Right].lastQRank, curQSymb};
            searchPart<OnMatchL, OnMatchR>(cursor, qcursorNext);
        }
        if (mismatchAllowed) {
            auto cursors = [&]() {
                if constexpr (Right) {
                    return cur.extendRight();
                } else {
                    return cur.extendLeft();
                }
            }();

            auto se = Restore{e, e+1};

            // search substitutes
            {
                auto sp = Restore{part, part+1};
                auto s1 = Restore{side[Right].lastRank};
                auto s2 = Restore{side[Right].lastQRank, curQSymb};
                for (size_t symb{1}; symb < Index::Sigma; ++symb) {
                    if (symb == curQSymb) continue;
                    side[Right].lastRank = symb;
                    searchPart<OnSubstituteL, OnSubstituteR>(cursors[symb],  qcursorNext);

                }
            }
            // search deletion
            if constexpr (Deletion) {
                auto s1 = Restore{side[Right].lastRank};

                for (size_t symb{1}; symb < Index::Sigma; ++symb) {
                    if (TInfo != 'M' || side[Right].lastQRank != symb) {
                        side[Right].lastRank = symb;
                        searchPart<OnDeletionL, OnDeletionR>(cursors[symb], qcur);
                    }

                }
            }

            // search insertion
            if constexpr (Insertion) {
                if (TInfo != 'M' || side[Right].lastQRank != curQSymb) {
                    auto sp = Restore{part, part+1};
                    auto s2 = Restore{side[Right].lastQRank, curQSymb};
                    searchPart<OnInsertionL, OnInsertionR>(cur, qcursorNext);
                }
            }
        }
    }
    #endif

    template <char LInfo, char RInfo, bool Right>
    void searchPartDirSingleIndex(cursor_t const& cur, qcursor_t const& qcur) const noexcept {
        static constexpr char TInfo = Right ? RInfo : LInfo;

        constexpr bool Deletion     = (TInfo != 'S' && TInfo != 'I') && Edit;
        constexpr bool Insertion    = (TInfo != 'S' && TInfo != 'D') && Edit;

        constexpr char OnMatchL      = Right ? LInfo : 'M';
        constexpr char OnMatchR      = Right ? 'M'   : RInfo;
        constexpr char OnSubstituteL = Right ? LInfo : 'S';
        constexpr char OnSubstituteR = Right ? 'S'   : RInfo;
        constexpr char OnDeletionL   = Right ? LInfo : 'D';
        constexpr char OnDeletionR   = Right ? 'D'   : RInfo;
        constexpr char OnInsertionL  = Right ? LInfo : 'I';
        constexpr char OnInsertionR  = Right ? 'I'   : RInfo;


        bool mismatchAllowed = l[part] <= e+1 and e+1 <= u[part];

        auto [curISymb, icursorNext] = [&]() -> std::tuple<size_t, cursor_t> {
            if constexpr (Right) {
                auto symb = cur.symbolRight();
                auto cur_  = cur.extendRight(symb);
                return {symb, cur_};
            } else {
                auto symb = cur.symbolLeft();
                auto cur_  = cur.extendLeft(symb);
                return {symb, cur_};
            }
        }();

        bool matchAllowed    = l[part] <= e and e <= u[part]
                               and ((TInfo != 'I' or curISymb != side[Right].lastQRank)
                               and (TInfo != 'D' or curISymb != side[Right].lastRank));

        if (mismatchAllowed) {
            auto qcursors = qcur.extendLeft();

            // search matches
            if (matchAllowed && curISymb != 0) {
                auto sp = Restore{part, part+1};
                auto s1 = Restore{side[Right].lastRank, curISymb};
                auto s2 = Restore{side[Right].lastQRank, curISymb};
                searchPart<OnMatchL, OnMatchR>(icursorNext, qcursors[curISymb]);
            }

            auto se = Restore{e, e+1};
            // search substitutes
            if (curISymb != 0) {
                auto s1 = Restore{side[Right].lastRank};
                auto s2 = Restore{side[Right].lastQRank};
                auto sp = Restore{part, part+1};
                for (size_t qsymb{1}; qsymb < QIndex::Sigma; ++qsymb) {
                    if (qcursors[qsymb].count() == 0) continue;
                    if (curISymb == qsymb) continue;
                    side[Right].lastQRank = qsymb;
                    side[Right].lastRank = curISymb;
                    searchPart<OnSubstituteL, OnSubstituteR>(icursorNext, qcursors[qsymb]);
                }
            }
            // search deletion
            if constexpr (Deletion) {
                if (curISymb != 0) {
                    auto s1 = Restore{side[Right].lastRank};
                    if (TInfo != 'M' || side[Right].lastQRank != curISymb) {
                        side[Right].lastRank = curISymb;
                        searchPart<OnDeletionL, OnDeletionR>(icursorNext, qcur);
                    }
                }
            }

            // search insertion
            if constexpr (Insertion) {
                auto s1 = Restore{side[Right].lastQRank};
                auto sp = Restore{part, part+1};
                auto lastQRank = side[Right].lastQRank;
                for (size_t qsymb{1}; qsymb < QIndex::Sigma; ++qsymb) {
                    if (TInfo != 'M' || lastQRank != qsymb) {
                        side[Right].lastQRank = qsymb;
                        searchPart<OnInsertionL, OnInsertionR>(cur, qcursors[qsymb]);
                    }
                }
            }
        } else if (matchAllowed && curISymb != 0) {
            auto qcursorNext = qcur.extendLeft(curISymb);
            auto sp = Restore{part, part+1};
            auto s1 = Restore{side[Right].lastRank, curISymb};
            auto s2 = Restore{side[Right].lastQRank, curISymb};
            searchPart<OnMatchL, OnMatchR>(icursorNext, qcursorNext);
        }
    }

};

template <bool Edit, typename Index, typename QIndex, typename search_t, typename delegate_t>
void search(Index const& index, QIndex const& queryIndex, search_t const& search, size_t threshold, size_t optMode, delegate_t&& delegate) {
    Search<Edit, Index, QIndex, search_t, delegate_t>{index, queryIndex, search, threshold, optMode, delegate};
}

}
