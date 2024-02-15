// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "SelectCursor.h"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

/**
 * As ng14
 * but direction change is predetermined
 */
namespace fmindex_collection::search_ng15 {

struct Block {
    size_t pi;
    size_t l;
    size_t u;
};

template <typename cursor_t, typename search_scheme_t, typename query_t, typename delegate_t, bool Right>
struct Search {
    constexpr static size_t Sigma = cursor_t::Sigma;

    using BlockIter = typename search_scheme_t::const_iterator;

    search_scheme_t const& search;
    query_t const& query;
    delegate_t const& delegate;

    Search(cursor_t const& _cursor, search_scheme_t const& _search, query_t const& _query, size_t e, delegate_t const& _delegate)
        : search    {_search}
        , query     {_query}
        , delegate  {_delegate}
    {
        auto blockIter = search.begin();

        search_next<'M'>(_cursor, e, blockIter);
    }

    static auto extend(cursor_t const& cur, uint8_t symb) noexcept {
        if constexpr (Right) {
            return cur.extendRight(symb);
        } else {
            return cur.extendLeft(symb);
        }
    }
    static auto extend(cursor_t const& cur) noexcept {
        if constexpr (Right) {
            return cur.extendRight();
        } else {
            return cur.extendLeft();
        }
    }

    template <char Info>
    void search_next(cursor_t const& cur, size_t e, BlockIter blockIter) const noexcept {

        if (cur.count() == 0) {
            return;
        }

        if (blockIter == end(search)) {
            if constexpr (Info == 'M' or Info == 'S') {
                delegate(cur, e);
            }
            return;
        }


        constexpr bool Deletion     = Info == 'M' or Info == 'D';
        constexpr bool Insertion    = Info == 'M' or Info == 'I';

        bool matchAllowed    = blockIter->l <= e and e <= blockIter->u;
        bool mismatchAllowed = blockIter->l <= e+1 and e+1 <= blockIter->u;

        auto symb = query[blockIter->pi];

        if (mismatchAllowed) {
            auto cursors = extend(cur);

            if (matchAllowed) {
                auto newCur = cursors[symb];
                search_next<'M'>(newCur, e, blockIter+1);
            }

            for (size_t i{1}; i < symb; ++i) {
                auto newCur = cursors[i];

                if constexpr (Deletion) {
                    search_next<'D'>(newCur, e+1, blockIter); // deletion occurred in query
                }
                search_next<'S'>(newCur, e+1, blockIter+1); // as substitution
            }

            for (size_t i(symb+1); i < Sigma; ++i) {
                auto newCur = cursors[i];

                if constexpr (Deletion) {
                    search_next<'D'>(newCur, e+1, blockIter); // deletion occurred in query
                }
                search_next<'S'>(newCur, e+1, blockIter+1); // as substitution
            }

            if constexpr (Insertion) {
                search_next<'I'>(cur, e+1, blockIter+1); // insertion occurred in query

            }

        } else if (matchAllowed) {
            auto newCur = extend(cur, symb);
            search_next<'M'>(newCur, e, blockIter+1);
        }
    }
};



template <typename index_t, typename queries_t, typename search_schemes_t, typename delegate_t>
void search(index_t const & index, queries_t && queries, search_schemes_t const & search_scheme, delegate_t && delegate)
{
    using cursor_t = select_cursor_t<index_t>;
    static_assert(not cursor_t::Reversed, "reversed fmindex is not supported");

    if (search_scheme.empty()) return;
    size_t qidx;
    auto internal_delegate = [&delegate, &qidx] (auto const & it, size_t e) {
        delegate(qidx, it, e);
    };


    std::vector<std::vector<std::vector<Block>>> search_scheme2;
    for (auto s : search_scheme) {
        std::vector<std::vector<Block>> search2;
        int lastDir = 0;
        for (size_t i{0}; i < s.pi.size(); ++i) {
            auto dir = [&]() {
                if (i == 0) {
                    return s.pi[i] < s.pi[i+1]?1:-1;
                } else {
                    return s.pi[i-1] < s.pi[i]?1:-1;
                }
            }();
            if (lastDir == 0) {
                search2.emplace_back();
                if (dir == -1) {
                    search2.emplace_back();
                }
            } else if (lastDir != dir) {
                search2.emplace_back();
            }
            lastDir = dir;
            search2.back().emplace_back(Block{size_t{s.pi[i]}, size_t{s.l[i]}, size_t{s.u[i]}});
        }
        search_scheme2.emplace_back(std::move(search2));
    }

    using Cursor = BiFMIndexCursor<index_t>;


    using Callback = std::function<void(Cursor const&, size_t e)>;
    Callback sch;

    using query_t = std::decay_t<decltype(queries[0])>;
    using search_t = std::decay_t<decltype(search_scheme2[0][0])>;

    query_t const* query;
    std::vector<search_t> const* search;
    size_t depth{};
    sch = [&](Cursor const& cursor, size_t e) {
        if (depth == search->size()) {
            internal_delegate(cursor, e);
            return;
        }
        if (depth % 2 == 0) {
            Search<Cursor, search_t, query_t, Callback, true>{cursor, search->at(depth++), *query, e, sch};
        } else {
            Search<Cursor, search_t, query_t, Callback, false>{cursor, search->at(depth++), *query, e, sch};
        }
        depth -= 1;
    };


    for (size_t i{0}; i < queries.size(); ++i) {
        qidx = i;
        query = &queries[qidx];
        for (size_t j{0}; j < search_scheme.size(); ++j) {
            search = &search_scheme2[j];
            // call sch
            sch(Cursor{index}, 0);
        }
    }

}

}
