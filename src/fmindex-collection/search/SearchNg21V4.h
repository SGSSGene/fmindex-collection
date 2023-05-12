// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "../BiFMIndexCursor.h"

#include <array>
#include <cstddef>

/**
 * like search_ng21V3 but no template functions
 */
namespace fmindex_collection {
namespace search_ng21V4 {

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

    Search(index_t const& _index, search_scheme_t const& _search, query_t const& _query, delegate_t const& _delegate)
        : index     {_index}
        , search    {_search}
        , query     {_query}
        , delegate  {_delegate}
    {
        auto cur       = cursor_t{index};

        search_next_dir(cur, 0, 'M', 'M', true);
    }

    static auto extend(cursor_t const& cur, uint8_t symb, bool right) noexcept {
        if (right) {
            return cur.extendRight(symb);
        } else {
            return cur.extendLeft(symb);
        }
    }
    static auto extend(cursor_t const& cur, bool right) noexcept {
        if (right) {
            return cur.extendRight();
        } else {
            return cur.extendLeft();
        }
    }


    void search_next(cursor_t const& cur, size_t lastRank, char LInfo, char RInfo) noexcept {
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
            search_next_dir(cur, lastRank, LInfo, RInfo, true);
        } else {
            cur.prefetchLeft();
            search_next_dir(cur, lastRank, LInfo, RInfo, false);
        }
    }

    void search_next_dir(cursor_t const& cur, size_t lastRank, char LInfo, char RInfo, bool right) noexcept {
        static char TInfo = right ? RInfo : LInfo;

        bool Deletion     = TInfo == 'M' or TInfo == 'D';
        bool Insertion    = TInfo == 'M' or TInfo == 'I';

        char OnMatchL      = right ? LInfo : 'M';
        char OnMatchR      = right ? 'M'   : RInfo;
        char OnSubstituteL = right ? LInfo : 'S';
        char OnSubstituteR = right ? 'S'   : RInfo;
        char OnDeletionL   = right ? LInfo : 'D';
        char OnDeletionR   = right ? 'D'   : RInfo;
        char OnInsertionL  = right ? LInfo : 'I';
        char OnInsertionR  = right ? 'I'   : RInfo;

        auto symb = query[search.pi[pos]];

        bool matchAllowed    = search.l[pos] <= e and e <= search.u[pos]
                               and (TInfo != 'I' or symb != query[search.pi[pos-1]])
                               and (TInfo != 'D' or symb != lastRank);
        bool mismatchAllowed = search.l[pos] <= e+1 and e+1 <= search.u[pos];

        if (mismatchAllowed) {
            auto cursors = extend(cur, right);

            pos += 1;
            if (matchAllowed) {
                search_next(cursors[symb], symb, OnMatchL, OnMatchR);
            }

            e+=1;
            for (uint8_t s{1}; s < Sigma; ++s) {
                if (s == symb) continue;
                search_next(cursors[s], s, OnSubstituteL, OnSubstituteR); // as substitution
            }

            if (Insertion) {
                search_next(cur, lastRank, OnInsertionL, OnInsertionR); // insertion occurred in query
            }
            pos -= 1;

            if (Deletion) {
                for (uint8_t s{1}; s < Sigma; ++s) {
                    if (s == symb) continue;
                    search_next(cursors[s], s, OnDeletionL, OnDeletionR); // deletion occurred in query
                }
            }

            e -= 1;
        } else if (matchAllowed) {
            auto newCur = extend(cur, symb, right);
            pos += 1;
            search_next(newCur, symb, OnMatchL, OnMatchR);
            pos -= 1;
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
}
