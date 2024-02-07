// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "SelectCursor.h"

namespace fmindex_collection::search_one_error {

/* Search algorithm with explicit programmed search scheme
 */
template <typename index_t, Sequence query_t, typename delegate_t>
struct Search {
    using cursor_t = BiFMIndexCursor<index_t>;
    constexpr static size_t Sigma = index_t::Sigma;

    index_t const& index;
    query_t const& query;
    delegate_t const& delegate;

    void search() {
        search_left_to_right();
        search_right_to_left();
    }

    void search_left_to_right() {
        // first half: left -> right
        auto cur = cursor_t{index};
        size_t mid = query.size() - query.size()/2;
        for (size_t i{0}; i < mid; ++i) {
            auto r = query[i];
            cur = cur.extendRight(r);
            if (cur.empty()) {
                return;
            }
        }
        search_left_to_right_second_half(cur, mid);
    }

    void search_left_to_right_second_half(cursor_t cur, size_t i) {
        // search second half: left -> right
        for (; i < query.size(); ++i) {
            auto r = query[i];
            auto nextCur = cur.extendRight();
            for (size_t s{1}; s < Sigma; ++s) {
                if (r != s) {
                    search_left_to_right_no_errors(nextCur[s], i+1);
                }
            }
            cur = nextCur[r];
            if (cur.empty()) {
                return;
            }
        }
        delegate(cur, 0);
    }

    void search_left_to_right_no_errors(cursor_t cur, size_t i) {
        // search second half: left -> right
        if (cur.empty()) return;

        for (; i < query.size(); ++i) {
            auto r = query[i];
            cur = cur.extendRight(r);
            if (cur.empty()) {
                return;
            }
        }
        delegate(cur, 1);
    }

    void search_right_to_left() {
        // first half: right -> left
        auto cur = cursor_t{index};
        size_t mid = query.size()/2;
        for (size_t i{0}; i < mid; ++i) {
            auto r = query[query.size() - i - 1];
            cur = cur.extendLeft(r);
            if (cur.empty()) {
                return;
            }
        }
        // searching second half
        search_right_to_left_second_half(cur, mid);
    }


    void search_right_to_left_second_half(cursor_t cur, size_t i) {
        // second half: right -> left
        for (; i < query.size(); ++i) {
            auto r = query[query.size() - i - 1];
            auto nextCur = cur.extendLeft();
            for (size_t s{1}; s < Sigma; ++s) {
                if (r != s) {
                    search_right_to_left_no_errors(nextCur[s], i+1);
                }
            }
            cur = nextCur[r];
            if (cur.empty()) {
                return;
            }
        }
    }

    void search_right_to_left_no_errors(cursor_t cur, size_t i) {
        // second half: right -> left
        if (cur.empty()) return;

        for (; i < query.size(); ++i) {
            auto r = query[query.size() - i - 1];
            cur = cur.extendLeft(r);
            if (cur.empty()) {
                return;
            }
        }
        delegate(cur, 1);
    }
};

#if __clang__ // broken clang 15, lets hope they catch on
template <typename index_t, Sequence query_t, typename delegate_t>
Search(index_t const&, query_t const&, delegate_t const&) -> Search<index_t, query_t, delegate_t>;
#endif

template <typename index_t, Sequence query_t, typename delegate_t>
void search(index_t const& index, query_t const& query, delegate_t&& delegate) {
    using cursor_t = select_cursor_t<index_t>;
    static_assert(not cursor_t::Reversed, "reversed fmindex is not supported");

    auto u = Search{index, query, delegate};
    u.search();
}

template <typename index_t, Sequences queries_t, typename delegate_t>
void search(index_t const& index, queries_t const& queries, delegate_t&& delegate) {
    using cursor_t = select_cursor_t<index_t>;
    static_assert(not cursor_t::Reversed, "reversed fmindex is not supported");

    for (size_t qidx{0}; qidx < queries.size(); ++qidx) {
        auto u = Search{index, queries[qidx], [&](auto cursor, auto error) {
            delegate(qidx, cursor, error);
        }};
        u.search();
    }
}

}
