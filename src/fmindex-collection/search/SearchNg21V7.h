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
 * like search_ng21V6 (with abort flag)
 * but using two buffers
 */
namespace fmindex_collection {
namespace search_ng21V7 {

enum class Dir : uint8_t { Left, Right };
template <typename T>
struct Block {
    T      rank;
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
    std::vector<std::vector<Block<size_t>>> searches;
    search_scheme_t const& search_scheme;
    delegate_t const& delegate;

    struct QueueEntry {
        std::vector<Block<size_t>> const& scheme;
        cursor_t cursor;
        size_t   pos;
        size_t   lastRank;

        using F = void (Search::*)(std::vector<Block<size_t>> const&, cursor_t const&, size_t e, size_t pos, size_t lastRank);
        F func;
    };

    struct Buffer {
        std::vector<QueueEntry> current{};
        std::vector<QueueEntry> after{};
    };

    Buffer buffer;

    using ReturnValue = std::decay_t<decltype(delegate(0, std::declval<cursor_t>(), 0))>;
    ReturnValue abort{};
    size_t ct{};
    size_t qidx{};


    Search(index_t const& _index, search_scheme_t const& _search_scheme, delegate_t const& _delegate)
        : index        {_index}
        , search_scheme{_search_scheme}
        , delegate     {_delegate}
    {

        // generate reordered searches
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
        searches = std::move(reordered);
    }

    template <typename query_t, typename bestHit_t>
    void search(size_t _qidx, query_t const& query, bestHit_t bestHit = false) {
        qidx = _qidx;
        abort = {};
        ct = 0;
        buffer.current.clear();
        buffer.after.clear();

        // initialize search schemes
        for (size_t j{0}; j < searches.size(); ++j) {
            auto& search = searches[j];
            for (size_t k {0}; k < search.size(); ++k) {
                search[k].rank = query[search_scheme[j].pi[k]];
            }
            search_next<'M', 'M'>(search, cursor_t{index}, 0, 0, 0);
        }

        // perform searches
        size_t e{1};
        while(!buffer.after.empty()) {
            std::swap(buffer.after, buffer.current);
            buffer.after.clear();
            for (auto const& q : buffer.current) {
                (this->*q.func)(q.scheme, q.cursor, e, q.pos, q.lastRank);
                if (abort) return;
            }
            if (bestHit and ct > 0) {
                break;
            }
            e += 1;
        }
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
    void search_next(std::vector<Block<size_t>> const& search, cursor_t const& cur, size_t e, size_t pos, size_t lastRank) {
        if (cur.count() == 0) {
            return;
        }

        if (pos == search.size()) {
            if constexpr ((LInfo == 'M' or LInfo == 'I') and (RInfo == 'M' or RInfo == 'I')) {
                ct += cur.count();
                abort = delegate(qidx, cur, e);
                return;
            }
            return;
        }
        if (search[pos].dir == Dir::Right) {
            cur.prefetchRight();
            search_next_dir<LInfo, RInfo, true>(search, cur, e, pos, lastRank);
        } else {
            cur.prefetchLeft();
            search_next_dir<LInfo, RInfo, false>(search, cur, e, pos, lastRank);
        }
    }

    template <char LInfo, char RInfo, bool Right>
    void search_next_dir(std::vector<Block<size_t>> const& search, cursor_t const& cur, size_t e, size_t pos, size_t lastRank) {
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

        auto blockIter = search.begin() + pos;
        auto symb = blockIter->rank;

        bool matchAllowed    = blockIter->l <= e and e <= blockIter->u
                               and (TInfo != 'I' or symb != (blockIter-1)->rank)
                               and (TInfo != 'D' or symb != lastRank);
        bool mismatchAllowed = blockIter->l <= e+1 and e+1 <= blockIter->u;

        if (mismatchAllowed) {
            auto cursors = extend<Right>(cur);

            if (matchAllowed and cursors[symb].count()) {
                search_next<OnMatchL, OnMatchR>(search, cursors[symb], e, pos+1, symb);
            }

#if __clang__
            for (uint8_t i{1}; i < symb; ++i) {
                if constexpr (Deletion) {
                    buffer.after.push_back(QueueEntry{search, cursors[i], pos, i, &Search::search_next<OnDeletionL, OnDeletionR>}); // deletion occurred in query
                }
                buffer.after.push_back(QueueEntry{search, cursors[i], pos+1, i, &Search::search_next<OnSubstituteL, OnSubstituteR>}); // as substitute
            }
            for (auto i{symb+1}; i < Sigma; ++i) {
                if constexpr (Deletion) {
                    buffer.after.push_back(QueueEntry{search, cursors[i], pos, i, &Search::search_next<OnDeletionL, OnDeletionR>}); // deletion occurred in query
                }
                buffer.after.push_back(QueueEntry{search, cursors[i], pos+1, i, &Search::search_next<OnSubstituteL, OnSubstituteR>}); // as substitute
            }

            if constexpr (Insertion) {
                buffer.after.push_back(QueueEntry{search, cur, pos+1, lastRank, &Search::search_next<OnInsertionL, OnInsertionR>}); // insertion occurred in query
            }

#else
            for (uint8_t i{1}; i < symb; ++i) {
                if constexpr (Deletion) {
                    buffer.after.emplace_back(search, cursors[i], pos, i, &Search::search_next<OnDeletionL, OnDeletionR>); // deletion occurred in query
                }
                buffer.after.emplace_back(search, cursors[i], pos+1, i, &Search::search_next<OnSubstituteL, OnSubstituteR>); // as substitute
            }
            for (auto i{symb+1}; i < Sigma; ++i) {
                if constexpr (Deletion) {
                    buffer.after.emplace_back(search, cursors[i], pos, i, &Search::search_next<OnDeletionL, OnDeletionR>); // deletion occurred in query
                }
                buffer.after.emplace_back(search, cursors[i], pos+1, i, &Search::search_next<OnSubstituteL, OnSubstituteR>); // as substitute
            }

            if constexpr (Insertion) {
                buffer.after.emplace_back(search, cur, pos+1, lastRank, &Search::search_next<OnInsertionL, OnInsertionR>); // insertion occurred in query
            }
#endif
        } else if (matchAllowed) {
            auto newCur = extend<Right>(cur, symb);
            if (newCur.count()) {
                search_next<OnMatchL, OnMatchR>(search, newCur, e, pos+1, symb);
            }
        }
    }
};


template <typename index_t, typename delegate_t>
auto refine_callback(delegate_t const& delegate) {
    using cursor_t = BiFMIndexCursor<index_t>;
    using R = std::decay_t<decltype(delegate(0, std::declval<cursor_t>(), 0))>;

    return [&](size_t qidx, auto cur, size_t e) {
        if constexpr (std::same_as<R, bool>) {
            return delegate(qidx, cur, e);
        } else {
            delegate(qidx, cur, e);
            return std::false_type{};
        }
    };
}

template <typename index_t, typename queries_t, typename search_scheme_t, typename delegate_t, typename bestHit_t = std::false_type>
void search(index_t const & index, queries_t && queries, search_scheme_t const & search_scheme, delegate_t && delegate, bestHit_t bestHit = {}) {
    using cursor_t = select_cursor_t<index_t>;
    static_assert(not cursor_t::Reversed, "reversed fmindex is not supported");

    auto internal_delegate = refine_callback<index_t>(delegate);

    auto search = Search{index, search_scheme, internal_delegate};
    for (size_t qidx{}; qidx < queries.size(); ++qidx) {
        search.search(qidx, queries[qidx], bestHit);
    }
}


template <typename index_t, typename queries_t, typename search_scheme_t, typename delegate_t, typename bestHit_t = std::false_type>
void search_n(index_t const & index, queries_t && queries, search_scheme_t const & search_scheme, size_t n, delegate_t && delegate, bestHit_t bestHit = {}) {
    using cursor_t = select_cursor_t<index_t>;
    static_assert(not cursor_t::Reversed, "reversed fmindex is not supported");

    size_t ct;
    auto cb = [&](size_t qidx, auto cur, size_t e) {
        if (cur.count() + ct > n) {
            cur.len = n-ct;
        }
        ct += cur.count();
        delegate(qidx, cur, e);
        return ct == n;
    };
    auto internal_delegate = refine_callback<index_t>(cb);

    auto search = Search{index, search_scheme, internal_delegate};

    for (size_t qidx{}; qidx < queries.size(); ++qidx) {
        ct = 0;
        search.search(qidx, queries[qidx], bestHit);
    }
}


template <typename index_t, typename queries_t, typename search_scheme_t, typename delegate_t>
void search_best(index_t const & index, queries_t && queries, search_scheme_t const & search_scheme, delegate_t && delegate) {
    return search(index, queries, search_scheme, delegate, std::true_type{});
}

template <typename index_t, typename queries_t, typename search_scheme_t, typename delegate_t>
void search_best_n(index_t const & index, queries_t && queries, search_scheme_t const & search_scheme, size_t n, delegate_t && delegate) {
    return search_n(index, std::forward<queries_t>(queries), search_scheme, n, std::forward<delegate_t>(delegate), std::true_type{});
}

}
}
