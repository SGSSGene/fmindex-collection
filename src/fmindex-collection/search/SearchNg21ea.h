// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "../BiFMIndexCursor.h"

#include <array>
#include <cstddef>

/**
 * Like Ng21, but aborts early
 * !TODO intext verification, configurable early abort
 */
namespace fmindex_collection {
namespace search_ng21ea {

enum class Dir : uint8_t { Left, Right };
template <typename T>
struct Block {
    T       rank;
    size_t l;
    size_t u;
    Dir dir;
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

        /*auto cursors = std::array<cursor_t, Sigma>{};
        for (size_t i{1}; i < Sigma; ++i) {
            if constexpr (Direction == Dir::Right) {
                cursors[i] = cur.extendRight(i);
            } else {
                cursors[i] = cur.extendLeft(i);
            }
        }
        return cursors;*/

    }


    template <char LInfo, char RInfo>
    void search_next(cursor_t const& cur, size_t e, BlockIter blockIter, size_t lastRank) const noexcept {
        if (cur.count() == 0) {
            return;
        }
        if (cur.count()*2 < size_t(std::distance(blockIter, end(search)))) {
            delegate(qidx, cur, e);
            return;
        }

        if (blockIter == end(search)) {
            if constexpr ((LInfo == 'M' or LInfo == 'I') and (RInfo == 'M' or RInfo == 'I')) {
                delegate(qidx, cur, e);
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
    void search_next_dir(cursor_t const& cur, size_t e, BlockIter blockIter, size_t lastRank) const noexcept {
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
                search_next<OnMatchL, OnMatchR>(newCur, e, blockIter+1, symb);
            }

            for (uint8_t i{1}; i < symb; ++i) {
                auto newCur = cursors[i];

                if constexpr (Deletion) {
                    search_next<OnDeletionL, OnDeletionR>(newCur, e+1, blockIter, i); // deletion occurred in query
                }
                search_next<OnSubstituteL, OnSubstituteR>(newCur, e+1, blockIter+1, i); // as substitution
            }

            for (uint8_t i(symb+1); i < Sigma; ++i) {
                auto newCur = cursors[i];

                if constexpr (Deletion) {
                    search_next<OnDeletionL, OnDeletionR>(newCur, e+1, blockIter, i); // deletion occurred in query
                }
                search_next<OnSubstituteL, OnSubstituteR>(newCur, e+1, blockIter+1, i); // as substitution
            }


            if constexpr (Insertion) {
                search_next<OnInsertionL, OnInsertionR>(cur, e+1, blockIter+1, lastRank); // insertion occurred in query
            }
        } else if (matchAllowed) {
            auto newCur = extend<Right>(cur, symb);
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
    auto internal_delegate = [&delegate] (size_t qidx, auto const & it, size_t e) {
        delegate(qidx, it, e);
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
}
