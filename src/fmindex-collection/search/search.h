// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "SearchPseudo.h"
#include "SearchNg21V2.h"
#include "SearchNoErrors.h"
#include "../search_scheme/generator/h2.h"
#include "../search_scheme/expand.h"

namespace fmc {

template <bool EditDistance, typename index_t, Sequence query_t, typename delegate_t>
void search(index_t const& _index, query_t const& _query, size_t _errors, delegate_t&& _delegate) {
    using cursor_t = select_cursor_t<index_t>;
    static_assert(not cursor_t::Reversed, "reversed fmindex is not supported");

    static thread_local auto cache = std::tuple<size_t, size_t, search_scheme::Scheme, search_scheme::Scheme>{std::numeric_limits<size_t>::max(), 0, {}, {}};
    // check if last scheme has correct errors and length, other wise generate it
    auto& [error, length, ss, search_scheme] = cache;
    if (error != _errors) { // regenerate everything
        ss            = search_scheme::generator::h2(_errors+2, 0, _errors);
        length        = _query.size();
        search_scheme = search_scheme::expand(ss, length);
    } else if (length != _query.size()) {
        length        = _query.size();
        search_scheme = search_scheme::expand(ss, length);
    }

    if constexpr (!EditDistance) {
        search_pseudo::search<EditDistance>(_index, _query, search_scheme, _delegate);
    } else {
        search_ng21V2::search(_index, _query, search_scheme, _delegate);
    }
}

template <typename index_t, Sequence query_t, typename delegate_t>
void search(index_t const& _index, query_t const& _query, delegate_t&& _delegate) {
    auto cursor = search_no_errors::search(_index, _query);
    if (cursor.count()) {
        _delegate(cursor);
    }
}

}
