// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../concepts.h"
#include "SelectCursor.h"


namespace fmc::search_no_errors {

template <typename index_t, Sequence query_t>
auto search(index_t const & index, query_t const& query) {
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
void search(index_t const & index, queries_t const& queries, delegate_t && delegate, size_t const BatchSize = 32) {
    using cursor_t = select_left_cursor_t<index_t>;
    static_assert(not cursor_t::Reversed, "reversed fmindex is not supported");


    auto queryIndices = std::vector<std::tuple<size_t, cursor_t>>{};
    size_t lastEntry = 0;

    auto fillBatch = [&]() {
        while (queryIndices.size() < BatchSize && lastEntry < queries.size()) {
            queryIndices.emplace_back(lastEntry, cursor_t{index});
            ++lastEntry;
        }
    };


    auto doBatchJump = [&]() {
        for (auto& [qidx, cur] : queryIndices) {
            auto const& query = queries[qidx];
            auto sym = query[query.size() - cur.steps - 1];
            cur = cur.extendLeft(sym);
        }
        queryIndices.erase(std::remove_if(begin(queryIndices), end(queryIndices), [&](auto const& v) {
            auto const& [qidx, cur] = v;
            auto const& query = queries[qidx];
            if (cur.count() == 0) return true;
            if (query.size() == cur.steps) {
                delegate(qidx, cur);
                return true;
            }
            return false;
        }), end(queryIndices));
    };

    fillBatch();
    while (!queryIndices.empty()) {
        doBatchJump();
        fillBatch();
    }
}

}
