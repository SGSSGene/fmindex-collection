#pragma once

#include "../concepts.h"

#include "../BiFMIndexCursor.h"


namespace fmindex_collection::search_no_errors {

template <typename index_t, Sequence query_t, typename delegate_t>
void search(index_t const & index, query_t && query, delegate_t && delegate) {
    using cursor_t = LeftBiFMIndexCursor<index_t>;
    static_assert(not cursor_t::Reversed, "reversed fmindex is not supported");

    auto cur = cursor_t{index};
    for (size_t i{0}; i < query.size(); ++i) {
        auto r = query[query.size() - i - 1];
        cur = cur.extendLeft(r);
        if (cur.empty()) {
            return;
        }
    }
    delegate(cur);
}

template <typename index_t, Sequences queries_t, typename delegate_t>
void search(index_t const & index, queries_t && queries, delegate_t && delegate) {

    for (size_t qidx{0}; qidx < queries.size(); ++qidx) {
        auto const& query = queries[qidx];
        search(index, query, [&qidx, &delegate](auto cursor) {
            delegate(qidx, cursor);
        });
    }
}

}
