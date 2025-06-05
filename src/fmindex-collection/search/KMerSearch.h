// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../concepts.h"
#include "SelectCursor.h"

namespace fmindex_collection::kmersearch {

template <typename index_t, Sequence query_t>
auto search(index_t const& index, query_t const& query) -> std::optional<select_left_cursor_t<index_t>> {
    using cursor_t = select_left_cursor_t<index_t>;
    static_assert(not cursor_t::Reversed, "reversed fmindex is not supported");

    auto cur = cursor_t{index};
    auto smallestCur = cur;
    for (size_t i{0}; i < query.size(); ++i) {
        auto s = query[query.size() - i - 1];
        cur = cur.extendLeft(s).clipToKMer();
        if (cur.count() == 0) return std::nullopt;
        if (cur.count() < smallestCur.count()) {
            smallestCur = cur;
        }
    }
    if (smallestCur.count() == index.size()) {
        return std::nullopt;
    }
    return smallestCur;
}

template <typename index_t, Sequences queries_t, typename delegate_t>
void search(index_t const& index, queries_t const& queries, delegate_t const& delegate) {

    for (size_t qidx{0}; qidx < queries.size(); ++qidx) {
        auto const& query = queries[qidx];
        auto cursor = search(index, query);
        if (cursor) {
            delegate(qidx, *cursor);
        }
    }
}

}
