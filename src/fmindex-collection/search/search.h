// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "SearchNg24.h"
#include "CachedSearchScheme.h"

namespace fmc {

template <bool EditDistance, typename index_t, Sequence query_t, typename delegate_t>
void search(index_t const& _index, query_t const& _query, size_t _errors, delegate_t&& _delegate) {
    search_ng24::search<EditDistance>(_index, _query, _errors, std::forward<delegate_t>(_delegate));
}

template <typename index_t, Sequence query_t, typename delegate_t>
void search(index_t const& _index, query_t const& _query, delegate_t&& _delegate) {
    auto cursor = search_no_errors::search(_index, _query);
    if (cursor.count()) {
        _delegate(cursor);
    }
}

}
