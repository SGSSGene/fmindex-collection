#pragma once

#include "../BiFMIndexCursor.h"

namespace fmindex_collection {
namespace search_one_error {

/* Search algorithm with explicit programmed search scheme
 */
template <typename index_t, typename queries_t, typename delegate_t>
struct Search {
    using cursor_t = BiFMIndexCursor<index_t>;
    constexpr static size_t Sigma = index_t::Sigma;
    using query_t = std::vector<uint8_t>;

    index_t const& index;
    queries_t const& queries;
    delegate_t const& delegate;

    void search() {
        for (size_t qidx{0}; qidx < queries.size(); ++qidx) {
            search_left_to_right_first_half(qidx, queries[qidx]);
            search_right_to_left_first_half(qidx, queries[qidx]);

        }
    }

    void search_right_to_left_first_half(size_t qidx, query_t const query) {
        auto cur = cursor_t{index};
        size_t mid = query.size()/2;
        for (size_t i{0}; i < mid; ++i) {
            auto r = query[query.size() - i - 1];
            cur = cur.extendLeft(r);
            if (cur.empty()) {
                return;
            }
        }
        search_right_to_left_second_half(qidx, query, cur, mid);
    }


    void search_right_to_left_second_half(size_t qidx, query_t const query, cursor_t cur, size_t i) {
        for (; i < query.size(); ++i) {
            auto r = query[query.size() - i - 1];
            auto nextCur = cur.extendLeft();
            for (size_t s{1}; s < Sigma; ++s) {
                if (r != s) {
                    search_right_to_left_no_errors(qidx, nextCur[s], query, i+1);
                }
            }
            cur = nextCur[r];
            if (cur.empty()) {
                return;
            }
        }
        delegate(qidx, cur, 0);
    }

    void search_right_to_left_no_errors(size_t qidx, cursor_t cur, query_t const& query, size_t i) {
        if (cur.empty()) return;

        for (; i < query.size(); ++i) {
            auto r = query[query.size() - i - 1];
            cur = cur.extendLeft(r);
            if (cur.empty()) {
                return;
            }
        }
        delegate(qidx, cur, 1);
    }

    void search_left_to_right_first_half(size_t qidx, query_t const query) {
        auto cur = cursor_t{index};
        size_t mid = query.size() - query.size()/2;
        for (size_t i{0}; i < mid; ++i) {
            auto r = query[i];
            cur = cur.extendRight(r);
            if (cur.empty()) {
                return;
            }
        }
        search_left_to_right_second_half(qidx, query, cur, mid);
    }

    void search_left_to_right_second_half(size_t qidx, query_t const query, cursor_t cur, size_t i) {
        for (; i < query.size(); ++i) {
            auto r = query[i];
            auto nextCur = cur.extendRight();
            for (size_t s{1}; s < Sigma; ++s) {
                if (r != s) {
                    search_left_to_right_no_errors(qidx, nextCur[s], query, i+1);
                }
            }
            cur = nextCur[r];
            if (cur.empty()) {
                return;
            }
        }
        delegate(qidx, cur, 0);
    }

    void search_left_to_right_no_errors(size_t qidx, cursor_t cur, query_t const& query, size_t i) {
        if (cur.empty()) return;

        for (; i < query.size(); ++i) {
            auto r = query[i];
            cur = cur.extendRight(r);
            if (cur.empty()) {
                return;
            }
        }
        delegate(qidx, cur, 1);
    }

};

template <typename index_t, typename queries_t, typename delegate_t>
void search(index_t const& index, queries_t&& queries, delegate_t&& delegate) {
    auto u = Search{index, queries, delegate};
    u.search();
}

}
}
