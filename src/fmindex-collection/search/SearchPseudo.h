// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <search_schemes/Scheme.h>

#include "../concepts.h"
#include "../fmindex/BiFMIndexCursor.h"
#include "SelectCursor.h"

namespace fmindex_collection::search_pseudo {

template <bool EditDistance, typename index_t, typename search_scheme_t, Sequence query_t, typename delegate_t>
struct Search {
    constexpr static size_t Sigma = index_t::Sigma;

    using cursor_t = select_cursor_t<index_t>;

    index_t const& index;

    decltype(search_scheme_t::pi) const& pi;
    decltype(search_scheme_t::l) const& l;
    decltype(search_scheme_t::u) const& u;

    query_t const& query;
    delegate_t const& delegate;

    Search(index_t const& _index, search_scheme_t const& _search, query_t const& _query, delegate_t const& _delegate) noexcept
        : index {_index}
        , pi{_search.pi}
        , l{_search.l}
        , u{_search.u}
        , query{_query}
        , delegate  {_delegate}
    {
        auto cur       = cursor_t{index};

        if constexpr (EditDistance) {
            search_distance(cur, 0, 0);
        } else {
            search_hm(cur, 0, 0);
        }
    }

    auto extend(cursor_t const& cur, uint8_t symb, std::size_t pos) const noexcept {
        if (pos == 0 or pi[pos-1] < pi[pos]) {
            return cur.extendRight(symb);
        } else {
            return cur.extendLeft(symb);
        }
    }
    auto extend(cursor_t const& cur, std::size_t pos) const noexcept {
        if (pos == 0 or pi[pos-1] < pi[pos]) {
            return cur.extendRight();
        } else {
            return cur.extendLeft();
        }
    }

    void search_hm(cursor_t const& cur, size_t e, std::size_t pos) const noexcept {
        if (cur.count() == 0) {
            return;
        }

        if (pos == query.size()) {
            if (l[pos-1] <= e and e <= u[pos-1]) {
                delegate(cur, e);
            }
            return;
        }

        if (e > u[pos]) {
            return;
        }


        // expected next character
        auto rank = query[pi[pos]];

        // cursors extended by one character
        auto cursors = std::array<cursor_t, Sigma>{};
        if (e+1 <= u[pos]) {
            for (size_t i{1}; i < Sigma; ++i) {
                cursors[i] = extend(cur, i, pos);
            }
        } else {
            cursors[rank] = extend(cur, rank, pos);
        }

        // search matches
        if (l[pos] <= e) {
            search_hm(cursors[rank], e,  pos+1);
        }

        // search substitution
        if (l[pos] <= e+1 and e+1 <= u[pos]) {
            for (size_t i{1}; i < Sigma; ++i) {
                if (i == rank) continue;
                search_hm(cursors[i], e+1, pos+1);
            }
        }
    }

    void search_distance(cursor_t const& cur, size_t e, std::size_t pos) const noexcept {

        if (cur.count() == 0) {
            return;
        }

        if (pos == query.size()) {
            if (l[pos-1] <= e and e <= u[pos-1]) {
                delegate(cur, e);
            }
            return;
        }

        if (e > u[pos]) {
            return;
        }

        // expected next character
        auto rank = query[pi[pos]];

        // cursors extended by one character
        auto cursors = std::array<cursor_t, Sigma>{};
        if (e+1 <= u[pos]) {
            for (size_t i{1}; i < Sigma; ++i) {
                cursors = extend(cur, pos);
            }
        } else {
            cursors[rank] = extend(cur, rank, pos);
        }

        // search matches
        if (l[pos] <= e) {
            search_distance(cursors[rank], e,  pos+1);
        }

        // search substitution
        if (l[pos] <= e+1 and e+1 <= u[pos]) {
            for (size_t i{1}; i < Sigma; ++i) {
                if (i == rank) continue;
                search_distance(cursors[i], e+1, pos+1);
            }
        }

        //search deletions
        if (e+1 <= u[pos]) {
            for (size_t i{1}; i < Sigma; ++i) {
                search_distance(cursors[i], e+1, pos);
            }
        }

        // search insertion
        if (l[pos] <= e+1 and e+1 <= u[pos]) {
            search_distance(cur, e+1, pos+1);
        }
    }
};


template <bool EditDistance, typename index_t, Sequences queries_t, typename search_schemes_t, typename delegate_t>
void search(index_t const & index, queries_t && queries, search_schemes_t const & search_scheme, delegate_t && delegate)
{
    using cursor_t = select_cursor_t<index_t>;
    static_assert(not cursor_t::Reversed, "reversed fmindex is not supported");

    std::size_t qidx;
    auto internal_delegate = [&qidx, &delegate] (auto const & it, size_t e) {
        delegate(qidx, it, e);
    };

    for (qidx = {0}; qidx < queries.size(); ++qidx) {
        for (size_t j{0}; j < search_scheme.size(); ++j) {
            Search<EditDistance, std::decay_t<decltype(index)>, std::decay_t<decltype(search_scheme[j])>, std::decay_t<decltype(queries[qidx])>, std::decay_t<decltype(internal_delegate)>> {index, search_scheme[j], queries[qidx], internal_delegate};
        }
    }
}
template <bool EditDistance, typename index_t, Sequence query_t, typename search_schemes_t, typename delegate_t>
void search(index_t const & index, query_t && query, search_schemes_t const & search_scheme, delegate_t && delegate)
{
    using cursor_t = select_cursor_t<index_t>;
    static_assert(not cursor_t::Reversed, "reversed fmindex is not supported");

    auto internal_delegate = [&delegate] (auto const & it, size_t e) {
        delegate(it, e);
    };

    for (size_t j{0}; j < search_scheme.size(); ++j) {
        Search<EditDistance, std::decay_t<decltype(index)>, std::decay_t<decltype(search_scheme[j])>, std::decay_t<decltype(query)>, std::decay_t<decltype(internal_delegate)>> {index, search_scheme[j], query, internal_delegate};
    }
}


}
