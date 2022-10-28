#pragma once

#include "SelectCursor.h"
#include "../concepts.h"

namespace fmindex_collection {
namespace search_backtracking_with_buffers {

/* Search algorithm with explicit programmed search scheme
 */
template <typename index_t, Sequence query_t, typename delegate_t>
struct Search {
    using cursor_t = select_cursor_t<index_t>;
    constexpr static size_t Sigma = index_t::Sigma;

    index_t const& index;
    query_t const& query;
    delegate_t const& delegate;

    size_t maxErrors;

    using buffer_t = std::vector<std::pair<cursor_t, size_t>>;
    buffer_t& buffer1;
    buffer_t& buffer2;

    void search() {
        auto cur = cursor_t{index};
        search_with_errors(cur);
    }

    template <typename ...Args>
    static auto extend(cursor_t cur, Args&&... args) {
        if constexpr (not cursor_t::Reversed) {
            return cur.extendLeft(std::forward<Args>(args)...);
        } else {
            return cur.extendRight(std::forward<Args>(args)...);
        }
    }

    void search_with_errors(cursor_t cur) {
        if (maxErrors == 0) {
            search_no_errors(cur, 0);
            return;
        }

        buffer1.emplace_back(cur, 0);

        for (size_t i{0}; i < query.size(); ++i) {
            auto r = query[query.size() - i - 1];

            for (auto & [cur, errorCount] : buffer1) {
                auto nextCur = extend(cur);
                for (size_t s{1}; s < Sigma; ++s) {
                    if (nextCur[s].empty()) continue;
                    auto newErrorCount = errorCount + ((r!=s)?1ul:0ul);
                    if (newErrorCount < maxErrors) {
                        buffer2.emplace_back(nextCur[s], newErrorCount);
                    } else {
                        search_no_errors(nextCur[s], i + 1);
                    }
                }
            }
            buffer1.clear();
            std::swap(buffer1, buffer2);
        }

        // last call must call all the delegates
        for (auto & [cur, errorCount] : buffer1) {
            delegate(cur, errorCount);
        }
        buffer1.clear();
    }

    void search_no_errors(cursor_t cur, size_t i) {
        for (; i < query.size(); ++i) {
            auto r = query[query.size() - i - 1];
            cur = extend(cur, r);
            if (cur.empty()) {
                return;
            }
        }
        delegate(cur, maxErrors);
    }
};

template <typename index_t, Sequence query_t, typename buffer_t, typename delegate_t>
void search(index_t const& index, query_t const& query, size_t maxError, buffer_t& buffer1, buffer_t& buffer2, delegate_t&& delegate) {
    auto u = Search{index, query, delegate, maxError, buffer1, buffer2};
    u.search();
}

}
}
