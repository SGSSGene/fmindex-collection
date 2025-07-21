// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "Restore.h"
#include "SelectCursor.h"

#include <array>
#include <cstddef>

/**
 * like search_ng24 but uses a scoring matrix
 */
namespace fmc::search_ng24sm {

/** Scoring matrix
 * sets identity to 0
 * removes first column and row (no matching possible)
 */
template <size_t QuerySigma, size_t RefSigma = QuerySigma>
struct ScoringMatrix {
    // 0: no pairing
    // 1: ambiguous match
    // 2: mismatch match
    std::array<std::array<size_t, RefSigma>, QuerySigma> matrix;

    std::array<std::vector<size_t>, QuerySigma> ambiguousList;
    std::array<std::vector<size_t>, QuerySigma> mismatchList;
    size_t allowedAmbiguous{};

    ScoringMatrix(size_t _allowedAmbiguous)
        : allowedAmbiguous{_allowedAmbiguous}
    {
        for (size_t y{1}; y < RefSigma; ++y) {
            for (size_t x{1}; x < QuerySigma; ++x) {
                setMismatch(x, y);
            }
        }
    }

    // sets the pairing as an ambiguous match
    void setAmbiguous(size_t queryRank, size_t refRank) {
        std::erase(ambiguousList[queryRank], refRank);
        std::erase(mismatchList[queryRank], refRank);
        matrix[queryRank][refRank] = 0;
        if (queryRank == refRank) return;
        matrix[queryRank][refRank] = 1;
        ambiguousList[queryRank].push_back(refRank);
    }

    // sets the pairing as a mismatch
    void setMismatch(size_t queryRank, size_t refRank) {
        std::erase(ambiguousList[queryRank], refRank);
        std::erase(mismatchList[queryRank], refRank);
        matrix[queryRank][refRank] = 0;
        if (queryRank == refRank) return;
        matrix[queryRank][refRank] = 2;
        mismatchList[queryRank].push_back(refRank);
    }

};


template <bool Edit, typename index_t, typename query_t, typename search_t, typename delegate_t>
struct Search {
    constexpr static size_t Sigma = index_t::Sigma;

    using cursor_t = select_cursor_t<index_t>;

    using SM = ScoringMatrix<Sigma>;

    index_t const& index;
    query_t const& query;
    search_t const& search;
    SM const& scoringMatrix;
    delegate_t const& delegate;

    struct Side {
        uint8_t lastRank{};
        uint8_t lastQRank{};
    };
    mutable std::array<Side, 2> side;
    mutable size_t e{};
    mutable size_t part{};
    mutable size_t remainingAmbiguous{};
    mutable size_t remainingMismatches{};

    Search(index_t const& _index, query_t const& _query, search_t const& _search, SM const& _scoringMatrix, size_t _maxErrors, delegate_t const& _delegate)
        : index         {_index}
        , query         {_query}
        , search        {_search}
        , scoringMatrix {_scoringMatrix}
        , delegate      {_delegate}
    {
        remainingAmbiguous  = std::min(_maxErrors, scoringMatrix.allowedAmbiguous);
        remainingMismatches = _maxErrors - remainingAmbiguous;
    }

    bool run() {
        auto cur       = cursor_t{index};
        return search_next<'M', 'M'>(cur);
    }


    template <bool Right>
    static auto extend(cursor_t const& cur, uint64_t symb) noexcept {
        if constexpr (Right) {
            return cur.extendRight(symb);
        } else {
            return cur.extendLeft(symb);
        }
    }
    template <bool Right>
    static auto extend(cursor_t const& cur) noexcept {
        if constexpr (Right) {
            return cur.extendRight();
        } else {
            return cur.extendLeft();
        }
    }


    template <char LInfo, char RInfo>
    bool search_next(cursor_t const& cur) const {
        if (cur.count() == 0) {
            return false;
        }

        if (part == search.pi.size()) {
            if constexpr (!Edit || ((LInfo == 'M' or LInfo == 'I') and (RInfo == 'M' or RInfo == 'I'))) {
                if (search.l[part-1] <= e and e <= search.u[part-1]) {
                    return delegate(cur, e);
                }
            }
            return false;
        }

        bool Right = (part==0) || (search.pi[part-1] < search.pi[part]);
        if (cur.count() > 1) {
            if (Right) {
                return search_next_dir<LInfo, RInfo, true>(cur);
            } else {
                return search_next_dir<LInfo, RInfo, false>(cur);
            }
        } else {
            if (Right) {
                return search_next_dir_single<LInfo, RInfo, true>(cur);
            } else {
                return search_next_dir_single<LInfo, RInfo, false>(cur);
            }
        }
    }

    template <char LInfo, char RInfo, bool Right>
    bool search_next_dir(cursor_t const& cur) const {
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

        auto nextSymb = query[search.pi[part]];

        bool matchAllowed    = search.l[part] <= e and e <= search.u[part]
                               and (TInfo != 'I' or nextSymb != side[Right].lastQRank)
                               and (TInfo != 'D' or nextSymb != side[Right].lastRank);
        bool mismatchAllowed = search.l[part] <= e+1 and e+1 <= search.u[part];

        if (mismatchAllowed) {
            auto cursors = extend<Right>(cur);

            if (matchAllowed) {
                auto const& newCur = cursors[nextSymb];
                auto s1 = Restore{side[Right].lastRank, nextSymb};
                auto s2 = Restore{side[Right].lastQRank, nextSymb};
                auto sp = Restore{part, part+1};
                bool f = search_next<OnMatchL, OnMatchR>(newCur);
                if (f) return true;
            }

            auto se = Restore{e, e+1};

            // search for as ambiguous substitutes
            if (remainingAmbiguous > 0) {
                auto s1 = Restore{side[Right].lastRank};
                auto s2 = Restore{side[Right].lastQRank, nextSymb};
                auto sp = Restore{part, part+1};
                auto sr = Restore{remainingAmbiguous, remainingAmbiguous-1};

                for (auto symb : scoringMatrix.ambiguousList[nextSymb]) {
                    auto const& newCur = cursors[symb];
                    side[Right].lastRank = symb;
                    bool f = search_next<OnSubstituteL, OnSubstituteR>(newCur); // as substitution
                    if (f) return true;
                }
            }

            if (remainingMismatches > 0) {
                auto sr = Restore{remainingMismatches, remainingMismatches-1};

                if constexpr (Deletion) {
                    auto s1 = Restore{side[Right].lastRank};
                    for (size_t symb{1}; symb < Sigma; ++symb) {
                        if (TInfo != 'M' || side[Right].lastQRank != symb) {
                            auto const& newCur = cursors[symb];
                            side[Right].lastRank = symb;
                            bool f = search_next<OnDeletionL, OnDeletionR>(newCur); // deletion occurred in query
                            if (f) return true;
                        }
                    }
                }
                {
                    auto s1 = Restore{side[Right].lastRank};
                    auto s2 = Restore{side[Right].lastQRank, nextSymb};
                    auto sp = Restore{part, part+1};

                    if (remainingAmbiguous == 0) { // if no ambiguous remaining, also search for them
                        for (auto symb : scoringMatrix.ambiguousList[nextSymb]) {
                            auto const& newCur = cursors[symb];
                            side[Right].lastRank = symb;
                            bool f = search_next<OnSubstituteL, OnSubstituteR>(newCur); // as substitution
                            if (f) return true;
                        }
                    }

                    for (auto symb : scoringMatrix.mismatchList[nextSymb]) {
                        auto const& newCur = cursors[symb];
                        side[Right].lastRank = symb;
                        bool f = search_next<OnSubstituteL, OnSubstituteR>(newCur); // as substitution
                        if (f) return true;
                    }
                }
            }

            if constexpr (Insertion) {
                if (remainingMismatches > 0) {
                    auto sr = Restore{remainingMismatches, remainingMismatches-1};
                    if (TInfo != 'M' || side[Right].lastQRank != nextSymb) {
                        auto s2 = Restore{side[Right].lastQRank, nextSymb};
                        auto sp = Restore{part, part+1};
                        bool f = search_next<OnInsertionL, OnInsertionR>(cur); // insertion occurred in query
                        if (f) return true;
                    }
                }
            }
        } else if (matchAllowed) {
            auto s1 = Restore{side[Right].lastRank, nextSymb};
            auto s2 = Restore{side[Right].lastQRank, nextSymb};
            auto newCur = extend<Right>(cur, nextSymb);
            auto sp = Restore{part, part+1};
            bool f = search_next<OnMatchL, OnMatchR>(newCur);
            if (f) return true;
        }
        return false;
    }

    template <char LInfo, char RInfo, bool Right>
    bool search_next_dir_single(cursor_t const& cur) const {
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


        auto [curISymb, icursorNext] = [&]() -> std::tuple<size_t, cursor_t> {
            if constexpr (Right) {
                auto symb = index.bwtRev.symbol(cur.lbRev);
                auto cur_  = cur.extendRight(symb);
                return {symb, cur_};
            } else {
                auto symb = index.bwt.symbol(cur.lb);
                auto cur_  = cur.extendLeft(symb);
                return {symb, cur_};
            }
        }();

        auto curQSymb = query[search.pi[part]];

        bool mismatchAllowed = search.l[part] <= e+1 and e+1 <= search.u[part];


        if constexpr (Insertion) {
            if (mismatchAllowed) {
                if (TInfo != 'M' || side[Right].lastQRank != curQSymb) {
                    auto se = Restore{e, e+1};
                    auto s2 = Restore{side[Right].lastQRank, curQSymb};
                    auto sp = Restore{part, part+1};

                    bool f = search_next<OnInsertionL, OnInsertionR>(cur);
                    if (f) return true;
                }
            }
        }




        // only insertions are possible
        if (curISymb == 0) {
            return false;
        }


        bool matchAllowed    = search.l[part] <= e and e <= search.u[part]
                               and (TInfo != 'I' or curQSymb != side[Right].lastQRank)
                               and (TInfo != 'D' or curQSymb != side[Right].lastRank);

        if (curISymb == curQSymb) {
            if (matchAllowed) {
                auto s1 = Restore{side[Right].lastRank, curQSymb};
                auto s2 = Restore{side[Right].lastQRank, curQSymb};
                auto sp = Restore{part, part+1};
                bool f = search_next<OnMatchL, OnMatchR>(icursorNext);
                if (f) return true;
            }
            if constexpr (Deletion) {
                if (mismatchAllowed && remainingMismatches > 0) {
                    auto sr = Restore{remainingMismatches, remainingMismatches-1};
                    if (TInfo != 'M' || side[Right].lastQRank != curISymb) {
                        auto se = Restore{e, e+1};
                        auto s1 = Restore{side[Right].lastRank, curISymb};
                        bool f = search_next<OnDeletionL, OnDeletionR>(icursorNext);
                        if (f) return true;
                    }
                }
            }
        } else if (mismatchAllowed) {
            auto se = Restore{e, e+1};

            // search substitute
            {
                auto s1 = Restore{side[Right].lastRank, curISymb};
                auto s2 = Restore{side[Right].lastQRank, curQSymb};
                auto sp = Restore{part, part+1};

                auto matchType = scoringMatrix.matrix[curQSymb][curISymb];
                if (matchType == 1) { // ambiguous match
                    if (remainingAmbiguous > 0) {
                        auto sr = Restore{remainingAmbiguous, remainingAmbiguous-1};
                        bool f = search_next<OnSubstituteL, OnSubstituteR>(icursorNext);
                        if (f) return true;
                    } else if (remainingMismatches > 0) {
                        auto sr = Restore{remainingMismatches, remainingMismatches-1};
                        bool f = search_next<OnSubstituteL, OnSubstituteR>(icursorNext);
                        if (f) return true;
                    }
                } else if (matchType == 2) { // mismatch
                    if (remainingMismatches > 0) {
                        auto sr = Restore{remainingMismatches, remainingMismatches-1};
                        bool f = search_next<OnSubstituteL, OnSubstituteR>(icursorNext);
                        if (f) return true;
                    }
                }
            }

            if constexpr (Deletion) {
                if (remainingMismatches > 0) {
                    auto sr = Restore{remainingMismatches, remainingMismatches-1};
                    if (TInfo != 'M' || side[Right].lastQRank != curISymb) {
                        auto s1 = Restore{side[Right].lastRank, curISymb};
                        bool f = search_next<OnDeletionL, OnDeletionR>(icursorNext);
                        if (f) return true;
                    }
                }
            }
        }
        return false;
    }

};

template <bool Edit=true, typename index_t, Sequence query_t, typename search_scheme_t, typename delegate_t>
void search(index_t const& index, query_t const& query, search_scheme_t const& search_scheme, ScoringMatrix<index_t::Sigma> const& scoringMatrix, delegate_t&& delegate) {
    using cursor_t = select_cursor_t<index_t>;
    using R = std::decay_t<decltype(delegate(std::declval<cursor_t>(), 0))>;

    auto internal_delegate = [&]() {
        if constexpr (std::same_as<R, bool>) {
            return delegate;
        } else {
            return [=](auto cur, size_t e) {
                delegate(cur, e);
                return std::false_type{};
            };
        }
    }();
    size_t maxErrors{};
    for (auto const& search : search_scheme) {
        maxErrors = std::max(search.u.back(), maxErrors);
    }


    for (auto const& search : search_scheme) {
        bool f = Search<Edit, index_t, query_t, decltype(search), decltype(internal_delegate)>{index, query, search, scoringMatrix, maxErrors, internal_delegate}.run();
        if (f) {
            return;
        }
    }
}

template <bool Edit=true, typename index_t, Sequences queries_t, typename search_scheme_t, typename delegate_t>
void search(index_t const& index, queries_t&& queries, search_scheme_t const& search_scheme, ScoringMatrix<index_t::Sigma> const& scoringMatrix, delegate_t&& delegate) {
    if (search_scheme.empty()) return;

    for (size_t qidx{}; qidx < queries.size(); ++qidx) {
        search<Edit>(index, queries[qidx], search_scheme, scoringMatrix, [&](auto const& cur, size_t e) {
            delegate(qidx, cur, e);
        });
    }
}

template <bool Edit=true, typename index_t, Sequence query_t, typename search_scheme_t, typename delegate_t>
void search_n(index_t const& index, query_t const& query, search_scheme_t const& search_scheme, size_t n, ScoringMatrix<index_t::Sigma> const& scoringMatrix, delegate_t&& delegate) {
    if (search_scheme.empty()) return;

    size_t ct{};
    search<Edit>(index, query, search_scheme, scoringMatrix, [&] (auto cur, size_t e) {
        if (cur.count() + ct > n) {
            cur.len = n-ct;
        }
        ct += cur.count();
        delegate(cur, e);
        return ct == n;
    });
}


template <bool Edit=true, typename index_t, Sequences queries_t, typename search_scheme_t, typename delegate_t>
void search_n(index_t const& index, queries_t&& queries, search_scheme_t const& search_scheme, size_t n, ScoringMatrix<index_t::Sigma> const& scoringMatrix, delegate_t&& delegate) {
    if (search_scheme.empty()) return;

    for (size_t qidx{}; qidx < queries.size(); ++qidx) {
        size_t ct{};
        search<Edit>(index, queries[qidx], search_scheme, scoringMatrix, [&] (auto cur, size_t e) {
            if (cur.count() + ct > n) {
                cur.len = n-ct;
            }
            ct += cur.count();
            delegate(qidx, cur, e);
            return ct == n;
        });
    }
}
}
