#pragma once

namespace fmindex_collection {
namespace search_no_errors {

template <typename index_t, typename queries_t, typename delegate_t>
void search(index_t const & index, queries_t && queries, delegate_t && delegate) {

    using cursor_t = LeftBiFMIndexCursor<index_t>;

    for (size_t qidx{0}; qidx < queries.size(); ++qidx) {
        auto const& query = queries[qidx];
        [&]() {
            auto cur = cursor_t{index};
            for (size_t i{0}; i < query.size(); ++i) {
                auto r = query[query.size() - i - 1];
                cur = cur.extendLeft(r);
                if (cur.empty()) {
                    return;
                }
            }
            delegate(qidx, cur, 0);
        }();
    }
}

}
}
