// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "SelectCursor.h"

#include <array>
#include <cstddef>

/**
 * like search_ng21
 * but using a abort flag instead of return values
 */
namespace fmindex_collection::search_ng21V6 {

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
    delegate_t const& delegate;

    using ReturnValue = std::decay_t<decltype(delegate(std::declval<cursor_t>(), 0))>;
    ReturnValue abort{};

    Search(index_t const& _index, search_scheme_t const& _search,  delegate_t const& _delegate)
        : index     {_index}
        , search    {_search}
        , delegate  {_delegate}
    {
        auto cur       = cursor_t{index};
        auto blockIter = search.begin();

        search_next<'M', 'M'>(cur, 0, blockIter, 0);
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
    void search_next(cursor_t const& cur, size_t e, BlockIter blockIter, size_t lastRank) {
        if (cur.count() == 0) {
            return;
        }

        if (blockIter == end(search)) {
            if constexpr ((LInfo == 'M' or LInfo == 'I') and (RInfo == 'M' or RInfo == 'I')) {
                abort = delegate(cur, e);
                return;
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
    void search_next_dir(cursor_t const& cur, size_t e, BlockIter blockIter, size_t lastRank) {
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
                if (abort) return;
            }

            for (size_t i{1}; i < symb; ++i) {
                auto newCur = cursors[i];

                if constexpr (Deletion) {
                    search_next<OnDeletionL, OnDeletionR>(newCur, e+1, blockIter, i); // deletion occurred in query
                    if (abort) return;
                }
                search_next<OnSubstituteL, OnSubstituteR>(newCur, e+1, blockIter+1, i); // as substitution
                if (abort) return;
            }

            for (size_t i(symb+1); i < Sigma; ++i) {
                auto newCur = cursors[i];

                if constexpr (Deletion) {
                    search_next<OnDeletionL, OnDeletionR>(newCur, e+1, blockIter, i); // deletion occurred in query
                    if (abort) return;
                }
                search_next<OnSubstituteL, OnSubstituteR>(newCur, e+1, blockIter+1, i); // as substitution
                if (abort) return;
            }


            if constexpr (Insertion) {
                search_next<OnInsertionL, OnInsertionR>(cur, e+1, blockIter+1, lastRank); // insertion occurred in query
                if (abort) return;
            }
        } else if (matchAllowed) {
            auto newCur = extend<Right>(cur, symb);
            search_next<OnMatchL, OnMatchR>(newCur, e, blockIter+1, symb);
            if (abort) return;
        }
    }
};



template <typename index_t, typename query_t, typename search_scheme_t, typename search_scheme_reordered_t, typename delegate_t>
void search_reordered(index_t const& index, query_t&& query, search_scheme_t const& search_scheme, search_scheme_reordered_t& reordered, delegate_t&& delegate) {
    using cursor_t = BiFMIndexCursor<index_t>;
    using R = std::decay_t<decltype(delegate(std::declval<cursor_t>(), 0))>;

    auto internal_delegate = [&]() {
        if constexpr (std::same_as<R, bool>) {
            return delegate;
        } else {
            return [=](auto cur, size_t e) {
                delegate(cur, e);
                return std::false_type{};
            };
        }
    }();

    for (size_t j{0}; j < search_scheme.size(); ++j) {
        auto& search = reordered[j];
        for (size_t k {0}; k < search.size(); ++k) {
            search[k].rank = query[search_scheme[j].pi[k]];
        }
        auto abort = Search{index, search, internal_delegate}.abort;
        if (abort) {
            return;
        }
    }
}

template <typename search_scheme_t>
auto prepare_reorder(search_scheme_t const& search_scheme) {
    auto reordered = std::vector<std::vector<Block<size_t>>>{};
    auto search2 = std::vector<Block<size_t>>{};
    for (auto const& s : search_scheme) {
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
        reordered.emplace_back(search2);
        search2.clear();
    }
    return reordered;
}

template <typename index_t, typename queries_t, typename search_scheme_t, typename delegate_t>
void search(index_t const & index, queries_t && queries, search_scheme_t const & search_scheme, delegate_t && delegate) {
    if (search_scheme.empty()) return;

    auto reordered = prepare_reorder(search_scheme);

    for (size_t qidx{}; qidx < queries.size(); ++qidx) {
        search_reordered(index, queries[qidx], search_scheme, reordered, [&](auto const& cur, size_t e) {
            delegate(qidx, cur, e);
        });
    }
}


template <typename index_t, typename queries_t, typename search_scheme_t, typename delegate_t>
void search_n(index_t const & index, queries_t && queries, search_scheme_t const & search_scheme, size_t n, delegate_t && delegate) {
    if (search_scheme.empty()) return;

    auto reordered = prepare_reorder(search_scheme);

    for (size_t qidx{}; qidx < queries.size(); ++qidx) {
        size_t ct{};
        search_reordered(index, queries[qidx], search_scheme, reordered, [&] (auto cur, size_t e) {
            if (cur.count() + ct > n) {
                cur.len = n-ct;
            }
            ct += cur.count();
            delegate(qidx, cur, e);
            return ct == n;
        });
    }
}

template <typename index_t, typename queries_t, typename search_scheme_t, typename delegate_t>
void search_best(index_t const & index, queries_t && queries, std::vector<search_scheme_t> const & search_schemes, delegate_t && delegate) {
    if (search_schemes.empty()) return;

    auto reordered_list = std::vector<decltype(prepare_reorder(search_schemes[0]))>{};
    for (auto const& search_scheme : search_schemes) {
        reordered_list.emplace_back(prepare_reorder(search_scheme));
    }

    for (size_t qidx{}; qidx < queries.size(); ++qidx) {
        for (size_t i{0}; i < reordered_list.size(); ++i) {
            auto& reordered     = reordered_list[i];
            auto& search_scheme = search_schemes[i];
            size_t ct{};
            search_reordered(index, queries[qidx], search_scheme, reordered, [&] (auto const& cur, size_t e) {
                ct += cur.count();
                delegate(qidx, cur, e);
            });
            if (ct > 0) break;
        }
    }
}

template <typename index_t, typename queries_t, typename search_scheme_t, typename delegate_t>
void search_best_n(index_t const & index, queries_t && queries, std::vector<search_scheme_t> const & search_schemes, size_t n, delegate_t && delegate) {
    using cursor_t = select_cursor_t<index_t>;
    static_assert(not cursor_t::Reversed, "reversed fmindex is not supported");

    if (search_schemes.empty()) return;

    auto reordered_list = std::vector<decltype(prepare_reorder(search_schemes[0]))>{};
    for (auto const& search_scheme : search_schemes) {
        reordered_list.emplace_back(prepare_reorder(search_scheme));
    }

    for (size_t qidx{}; qidx < queries.size(); ++qidx) {
        for (size_t i{0}; i < reordered_list.size(); ++i) {
            auto& reordered     = reordered_list[i];
            auto& search_scheme = search_schemes[i];
            size_t ct{};
            search_reordered(index, queries[qidx], search_scheme, reordered, [&] (auto cur, size_t e) {
                if (cur.count() + ct > n) {
                    cur.len = n-ct;
                }
                ct += cur.count();
                delegate(qidx, cur, e);
                return ct == n;
            });
            if (ct > 0) break;
        }
    }
}

}
