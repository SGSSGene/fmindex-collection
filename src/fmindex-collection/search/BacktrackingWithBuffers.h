// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
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
                    auto newErrorCount = errorCount + ((r!=s)?1ull:0ull);
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
#if __clang__ // !TODO broken clang 15, lets hope they catch on
template <typename index_t, Sequence query_t, typename delegate_t>
Search(index_t const&, query_t const&, delegate_t const&, size_t, auto&, auto&) -> Search<index_t, query_t, delegate_t>;
#endif

template <typename index_t, Sequence query_t, typename buffer_t, typename delegate_t>
void search(index_t const& index, query_t const& query, size_t maxError, buffer_t& buffer1, buffer_t& buffer2, delegate_t&& delegate) {
    using cursor_t = select_cursor_t<index_t>;
    if constexpr (not cursor_t::Reversed) {
        Search{index, query, delegate, maxError, buffer1, buffer2}.search();
    } else {
        #if __clang__ and __clang_major__ < 16 //!TODO some workaround, this is very expensive
            auto _query = query;
            std::ranges::reverse(_query);
            Search{index, _query, delegate, maxError, buffer1, buffer2}.search();
        #else
            Search{index, query | std::views::reverse, delegate, maxError, buffer1, buffer2}.search();
        #endif
    }
}

}
}
