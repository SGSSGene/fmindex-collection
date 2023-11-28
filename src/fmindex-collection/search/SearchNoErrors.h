// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../concepts.h"
#include "SelectCursor.h"


namespace fmindex_collection::search_no_errors {

template <typename index_t, Sequence query_t>
auto search(index_t const & index, query_t && query) {
    using cursor_t = select_left_cursor_t<index_t>;
    static_assert(not cursor_t::Reversed, "reversed fmindex is not supported");

    auto cur = cursor_t{index};
    for (size_t i{0}; i < query.size(); ++i) {
        auto r = query[query.size() - i - 1];
        cur = cur.extendLeft(r);
        if (cur.empty()) {
            return cur;
        }
    }
    return cur;
}

template <typename index_t, Sequences queries_t, typename delegate_t>
void search(index_t const & index, queries_t && queries, delegate_t && delegate) {

    for (size_t qidx{0}; qidx < queries.size(); ++qidx) {
        auto const& query = queries[qidx];
        auto cursor = search(index, query);
        delegate(qidx, cursor);
    }
}

}
