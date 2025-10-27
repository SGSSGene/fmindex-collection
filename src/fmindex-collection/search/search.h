// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../locate.h"
#include "SearchNg24.h"
#include "SearchNoErrors.h"
#include "CachedSearchScheme.h"

namespace fmc {

template <bool EditDistance, typename index_t, Sequence query_t, typename delegate_t>
void search(index_t const& _index, query_t const& _query, size_t _errors, delegate_t&& _delegate) {
    if (_errors == 0) {
        auto cursor = search_no_errors::search(_index, _query);
        if (cursor.count()) {
            _delegate(cursor, size_t{0});
        }
    } else {
        search_ng24::search<EditDistance>(_index, _query, _errors, std::forward<delegate_t>(_delegate));
    }
}

template <bool EditDistance, typename index_t, Sequences queries_t, typename delegate_t>
void search(index_t const& _index, queries_t const& _queries, size_t _errors, delegate_t&& _delegate) {
    if (_errors == 0) {
        search_no_errors::search(_index, _queries, [&](size_t qidx, auto const& cursor) {
            _delegate(qidx, cursor, size_t{0});
        });
    } else {
        search_ng24::search<EditDistance>(_index, _queries, _errors, std::forward<delegate_t>(_delegate));
    }
}

template <bool EditDistance, typename index_t, Sequence query_t, typename delegate_t>
void search_n(index_t const& _index, query_t const& _query, size_t _errors, size_t _n, delegate_t&& _delegate) {
    search_ng24::search<EditDistance>(_index, _query, _errors, std::forward<delegate_t>(_delegate), _n);
}

template <bool EditDistance, typename index_t, Sequences queries_t, typename delegate_t>
void search_n(index_t const& _index, queries_t const& _queries, size_t _errors, size_t _n, delegate_t&& _delegate) {
    search_ng24::search<EditDistance>(_index, _queries, _errors, std::forward<delegate_t>(_delegate), _n);
}

template <typename index_t, Sequences queries_t, typename delegate_t>
struct Search {
    index_t const&        index;
    queries_t const&      queries;
    bool                  editDistance{true};
    size_t                errors{0};
    std::optional<size_t> maxResults{};
    delegate_t const&     reportFunc;
    void operator()() {
        auto report = [&](size_t qidx, auto const& cursor, size_t errors) {
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;
                reportFunc(qidx, sid, spos+offset, errors);
            }
        };
        if (maxResults) {
            if (editDistance) {
                search_n<true>(index, queries, errors, *maxResults, report);
            } else {
                search_n<false>(index, queries, errors, *maxResults, report);
            }
        } else {
            if (editDistance) {
                search<true>(index, queries, errors, report);
            } else {
                search<false>(index, queries, errors, report);
            }
        }
    }
};

}
