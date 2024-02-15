// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "SelectCursor.h"

#include <array>
#include <cstddef>

/**
 * like search_ng21V3 but no template functions
 */
namespace fmindex_collection::search_ng21V5 {

template <typename index_t, typename search_scheme_t, typename query_t, typename delegate_t>
struct Search {
    constexpr static size_t Sigma = index_t::Sigma;

    using cursor_t = BiFMIndexCursor<index_t>;

    index_t const& index;
    search_scheme_t const& search;
    query_t const& query;
    delegate_t const& delegate;

    size_t e{};
    size_t pos{};

    char LInfo = 'M';
    char RInfo = 'M';

    Search(index_t const& _index, search_scheme_t const& _search, query_t const& _query, delegate_t const& _delegate)
        : index     {_index}
        , search    {_search}
        , query     {_query}
        , delegate  {_delegate}
    {
        auto cur       = cursor_t{index};

        search_next_dir_right(cur, 0);
    }

    void search_next(cursor_t const& cur, size_t lastRank) noexcept {
        if (cur.count() == 0) {
            return;
        }

        if (pos == query.size()) {
            if ((LInfo == 'M' or LInfo == 'I') and (RInfo == 'M' or RInfo == 'I')) {
                delegate(cur, e);
            }
            return;
        }
        auto goToRight = search.pi[pos-1] < search.pi[pos];

        if (goToRight) {
            cur.prefetchRight();
            search_next_dir_right(cur, lastRank);
        } else {
            cur.prefetchLeft();
            search_next_dir_left(cur, lastRank);
        }
    }
    void search_next_dir_right(cursor_t const& cur, size_t lastRank) noexcept {
	if (LInfo != 'M' and LInfo != 'I') return;
        auto symb = query[search.pi[pos]];

        bool matchAllowed    = search.l[pos] <= e and e <= search.u[pos]
                               and (RInfo != 'I' or symb != query[search.pi[pos-1]])
                               and (RInfo != 'D' or symb != lastRank);
        bool mismatchAllowed = search.l[pos] <= e+1 and e+1 <= search.u[pos];


        auto rInfo = RInfo;

        if (mismatchAllowed) {
            auto cursors = cur.extendRight();

            pos += 1;
            if (matchAllowed) {
                RInfo = 'M';
                search_next(cursors[symb], symb);
            }

            e+=1;
            RInfo = 'S';
            for (size_t s{1}; s < Sigma; ++s) {
                if (s == symb) continue;
                search_next(cursors[s], s); // as substitution
            }

            if (rInfo == 'M' or rInfo == 'I') {
                RInfo = 'I';
                search_next(cur, lastRank); // insertion occurred in query
            }
            pos -= 1;

            if (rInfo == 'M' or rInfo == 'D') {
                RInfo = 'D';
                for (size_t s{1}; s < Sigma; ++s) {
                    if (s == symb) continue;
                    search_next(cursors[s], s); // deletion occurred in query
                }
            }

            e -= 1;
        } else if (matchAllowed) {
            auto newCur = cur.extendRight(symb);
            RInfo = 'M';
            pos += 1;
            search_next(newCur, symb);
            pos -= 1;
        }
        RInfo = rInfo;
    }


    void search_next_dir_left(cursor_t const& cur, size_t lastRank) noexcept {
	if (RInfo != 'M' and RInfo != 'I') return;

        auto symb = query[search.pi[pos]];

        bool matchAllowed    = search.l[pos] <= e and e <= search.u[pos]
                               and (LInfo != 'I' or symb != query[search.pi[pos-1]])
                               and (LInfo != 'D' or symb != lastRank);
        bool mismatchAllowed = search.l[pos] <= e+1 and e+1 <= search.u[pos];


        auto lInfo = LInfo;

        if (mismatchAllowed) {
            auto cursors = cur.extendLeft();

            pos += 1;
            if (matchAllowed) {
                LInfo = 'M';
                search_next(cursors[symb], symb);
            }

            e+=1;
            LInfo = 'S';
            for (size_t s{1}; s < Sigma; ++s) {
                if (s == symb) continue;
                search_next(cursors[s], s); // as substitution
            }

            if (lInfo == 'M' or lInfo == 'I') {
                LInfo = 'I';
                search_next(cur, lastRank); // insertion occurred in query
            }
            pos -= 1;

            if (lInfo == 'M' or lInfo == 'D') {
                LInfo = 'D';
                for (size_t s{1}; s < Sigma; ++s) {
                    if (s == symb) continue;
                    search_next(cursors[s], s); // deletion occurred in query
                }
            }

            e -= 1;
        } else if (matchAllowed) {
            auto newCur = cur.extendLeft(symb);
            LInfo = 'M';
            pos += 1;
            search_next(newCur, symb);
            pos -= 1;
        }
        LInfo = lInfo;
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
