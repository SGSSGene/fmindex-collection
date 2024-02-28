// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "SelectCursor.h"

#include <array>
#include <cstddef>

/**
 * like search_ng21 but without extra vector
 */
namespace fmindex_collection::search_ng21V2 {

template <typename index_t, typename search_scheme_t, typename query_t, typename delegate_t>
struct Search {
    constexpr static size_t Sigma = index_t::Sigma;

    using cursor_t = BiFMIndexCursor<index_t>;

    index_t const& index;
    search_scheme_t const& search;
    query_t const& query;
    delegate_t const& delegate;

    Search(index_t const& _index, search_scheme_t const& _search, query_t const& _query, delegate_t const& _delegate)
        : index     {_index}
        , search    {_search}
        , query     {_query}
        , delegate  {_delegate}
    {
        auto cur       = cursor_t{index};

        search_next_dir<'M', 'M', true>(cur, 0, 0, 0);
    }

    template <bool Right>
    static auto extend(cursor_t const& cur, uint8_t symb) noexcept {
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
    void search_next(cursor_t const& cur, size_t e, size_t i, size_t lastRank) const noexcept {
        if (cur.count() == 0) {
            return;
        }

        if (i == query.size()) {
            if constexpr ((LInfo == 'M' or LInfo == 'I') and (RInfo == 'M' or RInfo == 'I')) {
                delegate(cur, e);
            }
            return;
        }
        auto goToRight = search.pi[i-1] < search.pi[i];

        if (goToRight) {
            cur.prefetchRight();
            search_next_dir<LInfo, RInfo, true>(cur, e, i, lastRank);
        } else {
            cur.prefetchLeft();
            search_next_dir<LInfo, RInfo, false>(cur, e, i, lastRank);
        }
    }

    template <char LInfo, char RInfo, bool Right>
    void search_next_dir(cursor_t const& cur, size_t e, size_t i, size_t lastRank) const noexcept {
        static constexpr char TInfo = Right ? RInfo : LInfo;

        constexpr bool Deletion     = TInfo == 'M' or TInfo == 'D';
        constexpr bool Insertion    = TInfo == 'M' or TInfo == 'I';

        constexpr char OnMatchL      = Right ? LInfo : 'M';
        constexpr char OnMatchR      = Right ? 'M'   : RInfo;
        constexpr char OnSubstituteL = Right ? LInfo : 'S';
        constexpr char OnSubstituteR = Right ? 'S'   : RInfo;
        constexpr char OnDeletionL   = Right ? LInfo : 'D';
        constexpr char OnDeletionR   = Right ? 'D'   : RInfo;
        constexpr char OnInsertionL  = Right ? LInfo : 'I';
        constexpr char OnInsertionR  = Right ? 'I'   : RInfo;

        auto symb = query[search.pi[i]];

        bool matchAllowed    = search.l[i] <= e and e <= search.u[i]
                               and (TInfo != 'I' or symb != query[search.pi[i-1]])
                               and (TInfo != 'D' or symb != lastRank);
        bool mismatchAllowed = search.l[i] <= e+1 and e+1 <= search.u[i];

        if (mismatchAllowed) {
            auto cursors = extend<Right>(cur);

            if (matchAllowed) {
                auto newCur = cursors[symb];
                search_next<OnMatchL, OnMatchR>(newCur, e, i+1, symb);
            }

            for (size_t s{1}; s < symb; ++s) {
                auto newCur = cursors[s];

                if constexpr (Deletion) {
                    search_next<OnDeletionL, OnDeletionR>(newCur, e+1, i, s); // deletion occurred in query
                }
                search_next<OnSubstituteL, OnSubstituteR>(newCur, e+1, i+1, s); // as substitution
            }

            for (size_t s(symb+1); s < Sigma; ++s) {
                auto newCur = cursors[s];

                if constexpr (Deletion) {
                    search_next<OnDeletionL, OnDeletionR>(newCur, e+1, i, s); // deletion occurred in query
                }
                search_next<OnSubstituteL, OnSubstituteR>(newCur, e+1, i+1, s); // as substitution
            }


            if constexpr (Insertion) {
                search_next<OnInsertionL, OnInsertionR>(cur, e+1, i+1, lastRank); // insertion occurred in query
            }
        } else if (matchAllowed) {
            auto newCur = extend<Right>(cur, symb);
            search_next<OnMatchL, OnMatchR>(newCur, e, i+1, symb);
        }
    }
};


template <typename index_t, typename queries_t, typename search_schemes_t, typename delegate_t>
void search(index_t const & index, queries_t && queries, search_schemes_t const & search_scheme, delegate_t && delegate)
{
    using cursor_t = select_cursor_t<index_t>;
    static_assert(not cursor_t::Reversed, "reversed fmindex is not supported");

    size_t qidx;
    auto internal_delegate = [&delegate, &qidx] (auto const & it, size_t e) {
        delegate(qidx, it, e);
    };

    for (qidx = 0; qidx < queries.size(); ++qidx) {
        auto const& query = queries[qidx];
        for (size_t j{0}; j < search_scheme.size(); ++j) {
            auto search = search_scheme[j];
            Search{index, search, query, internal_delegate};
        }
    }

}

}
