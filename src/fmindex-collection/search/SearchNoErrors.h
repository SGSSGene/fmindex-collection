#pragma once

namespace fmindex_collection {
namespace search_no_errors {

template <typename index_t, typename queries_t, typename delegate_t>
void search(index_t const & index, queries_t && queries, delegate_t && delegate) {

    using cursor_t = BiFMIndexCursor<index_t>;

    for (size_t qidx{0}; qidx < queries.size(); ++qidx) {
        auto const& query = queries[qidx];
        [&]() {
            auto cur = cursor_t{index};
            for (auto r : query) {
                cur = cur.extendRight(r);
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
