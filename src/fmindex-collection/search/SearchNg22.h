// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "SelectCursor.h"

#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <vector>


/**
 * like search_ng21 but with alignment string
 */
namespace fmindex_collection::search_ng22 {

enum class Dir : uint8_t { Left, Right };
template <typename T>
struct Block {
    T       rank;
    size_t l;
    size_t u;
    Dir dir;
};

template <bool Right>
struct ActionScope {
    std::deque<char>& list;
    ActionScope() = delete;
    ActionScope(ActionScope const&) = delete;
    ActionScope(ActionScope&&) = delete;
    ActionScope(std::deque<char>& list, char action)
        : list{list}
    {
        if constexpr (Right) {
            list.emplace_back(action);
        } else {
            list.emplace_front(action);
        }
    }
    ~ActionScope() {
        if constexpr (Right) {
            list.pop_back();
        } else {
            list.pop_front();
        }
    }
    ActionScope& operator=(ActionScope const&) = delete;
    ActionScope& operator=(ActionScope&&) = delete;
};

struct ActionList {
    std::deque<char> list;

    template <bool Right>
    auto push(char c) {
        return ActionScope<Right>{list, c};
    }
};

template <typename index_t, typename search_scheme_t, typename delegate_t>
struct Search {
    constexpr static size_t Sigma = index_t::Sigma;

    using cursor_t = BiFMIndexCursor<index_t>;
    using BlockIter = typename search_scheme_t::const_iterator;

    index_t const& index;
    search_scheme_t const& search;
    size_t qidx;
    delegate_t const& delegate;
    ActionList actions;

    Search(index_t const& _index, search_scheme_t const& _search, size_t _qidx, delegate_t const& _delegate)
        : index     {_index}
        , search    {_search}
        , qidx      {_qidx}
        , delegate  {_delegate}
    {
        auto cur       = cursor_t{index};
        auto blockIter = search.begin();

        search_next<'M', 'M'>(cur, 0, blockIter,0 );
    }

    template <bool Right>
    static auto extend(cursor_t const& cur, uint8_t symb) noexcept {
        if constexpr (Right) {
            return cur.extendRight(symb);
        } else {
            return cur.extendLeft(symb);
        }
    }
    template <bool Right>
    static auto extend(cursor_t const& cur) noexcept {
        if constexpr (Right) {
            return cur.extendRight();
        } else {
            return cur.extendLeft();
        }
    }


    template <char LInfo, char RInfo>
    void search_next(cursor_t const& cur, size_t e, BlockIter blockIter, size_t lastRank) noexcept {
        if (cur.count() == 0) {
            return;
        }

        if (blockIter == end(search)) {
            if constexpr ((LInfo == 'M' or LInfo == 'S') and (RInfo == 'M' or RInfo == 'S')) {
                delegate(qidx, cur, e, actions.list);
            }
            return;
        }
        if (blockIter->dir == Dir::Right) {
            cur.prefetchRight();
            search_next_dir<LInfo, RInfo, true>(cur, e, blockIter, lastRank);
        } else {
            cur.prefetchLeft();
            search_next_dir<LInfo, RInfo, false>(cur, e, blockIter, lastRank);
        }
    }

    template <char LInfo, char RInfo, bool Right>
    void search_next_dir(cursor_t const& cur, size_t e, BlockIter blockIter, size_t lastRank) noexcept {
        static constexpr char TInfo = Right ? RInfo : LInfo;

        constexpr bool Deletion     = TInfo == 'M' or TInfo == 'D';
        constexpr bool Insertion    = TInfo == 'M' or TInfo == 'I';

        constexpr char OnMatchL      = Right ? LInfo : 'M';
        constexpr char OnMatchR      = Right ? 'M'   : RInfo;
        constexpr char OnSubstituteL = Right ? LInfo : 'S';
        constexpr char OnSubstituteR = Right ? 'S'   : RInfo;
        constexpr char OnDeletionL   = Right ? LInfo : 'D';
        constexpr char OnDeletionR   = Right ? 'D'   : RInfo;
        constexpr char OnInsertionL  = Right ? LInfo : 'I';
        constexpr char OnInsertionR  = Right ? 'I'   : RInfo;

        auto symb = blockIter->rank;

        bool matchAllowed    = blockIter->l <= e and e <= blockIter->u
                               and (TInfo != 'I' or symb != (blockIter-1)->rank)
                               and (TInfo != 'D' or symb != lastRank);
        bool mismatchAllowed = blockIter->l <= e+1 and e+1 <= blockIter->u;

        if (mismatchAllowed) {
            auto cursors = extend<Right>(cur);

            if (matchAllowed) {
                auto newCur = cursors[symb];
                auto s = actions.push<Right>('M');
                search_next<OnMatchL, OnMatchR>(newCur, e, blockIter+1, symb);
            }

            for (size_t i{1}; i < symb; ++i) {
                auto newCur = cursors[i];

                if constexpr (Deletion) {
                    auto s = actions.push<Right>('D');
                    search_next<OnDeletionL, OnDeletionR>(newCur, e+1, blockIter, i); // deletion occurred in query
                }
                auto s = actions.push<Right>('S');
                search_next<OnSubstituteL, OnSubstituteR>(newCur, e+1, blockIter+1, i); // as substitution
            }

            for (size_t i(symb+1); i < Sigma; ++i) {
                auto newCur = cursors[i];

                if constexpr (Deletion) {
                    auto s = actions.push<Right>('D');
                    search_next<OnDeletionL, OnDeletionR>(newCur, e+1, blockIter, i); // deletion occurred in query
                }
                auto s = actions.push<Right>('S');
                search_next<OnSubstituteL, OnSubstituteR>(newCur, e+1, blockIter+1, i); // as substitution
            }

            if constexpr (Insertion) {
                auto s = actions.push<Right>('I');
                search_next<OnInsertionL, OnInsertionR>(cur, e+1, blockIter+1, lastRank); // insertion occurred in query
            }
        } else if (matchAllowed) {
            auto newCur = extend<Right>(cur, symb);
            auto s = actions.push<Right>('M');
            search_next<OnMatchL, OnMatchR>(newCur, e, blockIter+1, symb);
        }
    }
};


template <typename index_t, typename queries_t, typename search_schemes_t, typename delegate_t>
void search(index_t const & index, queries_t && queries, search_schemes_t const & search_scheme, delegate_t && delegate)
{
    using cursor_t = select_cursor_t<index_t>;
    static_assert(not cursor_t::Reversed, "reversed fmindex is not supported");

    if (search_scheme.empty()) return;
    auto internal_delegate = [&delegate] (size_t qidx, auto const & it, size_t e, auto const& actions) {
        delegate(qidx, it, e, actions);
    };

    std::vector<std::vector<Block<size_t>>> search_scheme2;
    for (auto s : search_scheme) {
        std::vector<Block<size_t>> search2;
        for (size_t i{0}; i < s.pi.size(); ++i) {
            auto dir = [&]() {
                if (i == 0) {
                    return s.pi[i] < s.pi[i+1]?Dir::Right:Dir::Left;
                } else {
                    return s.pi[i-1] < s.pi[i]?Dir::Right:Dir::Left;
                }
            }();
            search2.emplace_back(Block<size_t>{{}, size_t{s.l[i]}, size_t{s.u[i]}, dir});
        }
        search_scheme2.emplace_back(std::move(search2));
    }

    for (size_t i{0}; i < queries.size(); ++i) {
        auto const& query = queries[i];
        for (size_t j{0}; j < search_scheme.size(); ++j) {
            auto& search = search_scheme2[j];
            for (size_t k {0}; k < search.size(); ++k) {
                search[k].rank = query[search_scheme[j].pi[k]];
            }
            Search<std::decay_t<decltype(index)>, std::decay_t<decltype(search)>, std::decay_t<decltype(internal_delegate)>>{index, search, i, internal_delegate};
        }
    }

}

}
