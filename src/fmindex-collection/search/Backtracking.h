#pragma once

#include "SelectCursor.h"

namespace fmindex_collection {
namespace search_backtracking {

/* Search algorithm with explicit programmed search scheme
 */
template <typename index_t, typename queries_t, typename delegate_t>
struct Search {
    using cursor_t = select_cursor_t<index_t>;
    constexpr static size_t Sigma = index_t::Sigma;
    using query_t = std::vector<uint8_t>;

    index_t const& index;
    queries_t const& queries;
    delegate_t const& delegate;

    size_t maxErrors{};

    void search() {
        auto cur = cursor_t{index};
        for (size_t qidx{0}; qidx < queries.size(); ++qidx) {
            search_with_errors(0, qidx, queries[qidx], cur, 0);
        }
    }

    template <typename ...Args>
    static auto extend(cursor_t cur, Args&&... args) {
        if constexpr (not cursor_t::Reversed) {
            return cur.extendLeft(std::forward<Args>(args)...);
        } else {
            return cur.extendRight(std::forward<Args>(args)...);
        }
    }

    void search_with_errors(size_t e, size_t qidx, query_t const query, cursor_t cur, size_t i) {
        if (cur.empty()) {
            return;
        }
        if (e == maxErrors) {
            search_no_errors(qidx, query, cur, i);
            return;
        }
        for (; i < query.size(); ++i) {
            auto r = query[query.size() - i - 1];
            auto nextCur = extend(cur);
            for (size_t s{1}; s < Sigma; ++s) {
                if (r != s) {
                    search_with_errors(e+1, qidx, query, nextCur[s], i+1);
                }
            }
            cur = nextCur[r];
            if (cur.empty()) {
                return;
            }
        }
        delegate(qidx, cur, e);
    }

    void search_no_errors(size_t qidx, query_t const& query, cursor_t cur, size_t i) {
        if (cur.empty()) return;

        for (; i < query.size(); ++i) {
            auto r = query[query.size() - i - 1];
            cur = extend(cur, r);
            if (cur.empty()) {
                return;
            }
        }
        delegate(qidx, cur, maxErrors);
    }
};

template <typename index_t, typename queries_t, typename delegate_t>
void search(index_t const& index, queries_t&& queries, size_t maxError, delegate_t&& delegate) {
    auto u = Search{index, queries, delegate, maxError};
    u.search();
}

}
}
