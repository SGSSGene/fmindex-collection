// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "CachedSearchScheme.h"
#include "Restore.h"
#include "SelectCursor.h"

#include <array>
#include <cstddef>

/**
 * like search_ng25 but:
 *  - removes all the templates from the functions
 */
namespace fmc::search_ng26 {

template <bool Edit, typename index_t, typename query_t, typename search_t, typename delegate_t>
struct Search {
    constexpr static size_t Sigma = index_t::Sigma;
    constexpr static size_t FirstSymb = []() -> size_t {
        if constexpr (requires() { { index_t::FirstSymb }; }) {
            return index_t::FirstSymb;
        }
        return 1;
    }();

    using cursor_t = select_cursor_t<index_t>;

    index_t const& index;
    query_t const& query;
    search_t const& search;
    std::vector<size_t> const& partition;
    delegate_t const& delegate;

    struct Side {
        uint8_t lastRank{};
        uint8_t lastQRank{};
    };

    struct State {
        cursor_t cur;
        std::array<Side, 2> side;
        size_t e{};
        size_t part{};
        size_t partitionEntryValue;
        size_t queryPosL{};
        size_t queryPosR{};
        char LInfo, RInfo;
        bool Right{};
        bool NextPos{};
    };

    Search(index_t const& _index, query_t const& _query, search_t const& _search, std::vector<size_t> const& _partition, delegate_t const& _delegate)
        : index     {_index}
        , query     {_query}
        , search    {_search}
        , partition {_partition}
        , delegate  {_delegate}
    {}

    bool run() {
        auto state = State{};
        // check how many characters are before the first query char
        for (size_t i{0}; i < search.pi[0]; ++i) {
            state.queryPosL += partition[i];
            state.queryPosR += partition[i];
        }
        // Notice: queryPosL might 'underflow' (also at other spots)
        // but that is fine, after an underflow this variable is not
        // being accessed before overflowing again.
        state.queryPosL -= 1;
        state.partitionEntryValue = partition[search.pi[0]];
        state.cur = cursor_t{index};
        state.LInfo = 'M';
        state.RInfo = 'M';

        return search_next(state);
    }


    auto extend(State const& state, uint64_t symb) const noexcept {
        if (state.Right) {
            return state.cur.extendRight(symb);
        } else {
            return state.cur.extendLeft(symb);
        }
    }
    auto extend(State const& state) const noexcept {
        if (state.Right) {
            return state.cur.extendRight();
        } else {
            return state.cur.extendLeft();
        }
    }


    bool search_next(State const& state) const {
        if (state.cur.count() == 0) {
            return false;
        }
        if (state.part == partition.size()) {
            if (!Edit || ((state.LInfo == 'M' or state.LInfo == 'I') and (state.RInfo == 'M' or state.RInfo == 'I'))) {
                if (search.l.back() <= state.e and state.e <= search.u.back()) {
                    return delegate(state.cur, state.e);
                }
            }
            return false;
        }

        auto newState = state;
        newState.Right = (state.part==0) || (search.pi[state.part-1] < search.pi[state.part]);
        if (state.cur.count() > 1) {
            return search_next_dir(newState);
        } else {
            return search_next_dir_single(newState);
        }
    }

    bool search_next_pos(State state) const {
        if (state.cur.count() == 0) return false;
        if (state.NextPos) {
            if (state.Right) state.queryPosR += 1;
            else state.queryPosL -= 1;
            state.partitionEntryValue -= 1;

            if (state.partitionEntryValue == 0) {
                state.part += 1;
                if (state.part != partition.size()) {
                    state.partitionEntryValue = partition[search.pi[state.part]];
                }
                bool res = search_next(state);
                return res;
            }
        }

        if (state.cur.count() > 1) {
            return search_next_dir(state);
        } else {
            return search_next_dir_single(state);
        }
    }

    bool search_next_dir(State const& state) const {
        char const TInfo = state.Right ? state.RInfo : state.LInfo;

        bool const Deletion     = (TInfo != 'S' && TInfo != 'I') && Edit;
        bool const Insertion    = (TInfo != 'S' && TInfo != 'D') && Edit;

        char const OnMatchL      = state.Right ? state.LInfo : 'M';
        char const OnMatchR      = state.Right ? 'M'   : state.RInfo;
        char const OnSubstituteL = state.Right ? state.LInfo : 'S';
        char const OnSubstituteR = state.Right ? 'S'   : state.RInfo;
        char const OnDeletionL   = state.Right ? state.LInfo : 'D';
        char const OnDeletionR   = state.Right ? 'D'   : state.RInfo;
        char const OnInsertionL  = state.Right ? state.LInfo : 'I';
        char const OnInsertionR  = state.Right ? 'I'   : state.RInfo;

        auto nextSymb = query[state.Right?state.queryPosR:state.queryPosL];

        bool matchAllowed    = (state.partitionEntryValue > 1 or search.l[state.part] <= state.e)
                               and state.e <= search.u[state.part]
                               and (TInfo != 'I' or nextSymb != state.side[state.Right].lastQRank)
                               and (TInfo != 'D' or nextSymb != state.side[state.Right].lastRank);
        bool insertionAllowed    = (state.partitionEntryValue > 1 or search.l[state.part] <= state.e+1)
                                   and state.e+1 <= search.u[state.part];
        bool substitutionAllowed = insertionAllowed;
        bool mismatchAllowed     = state.e+1 <= search.u[state.part];

        if (mismatchAllowed) {
            auto cursors = extend(state);

            if (matchAllowed) {
                auto newState = state;
                newState.cur = cursors[nextSymb];
                newState.side[state.Right].lastRank = nextSymb;
                newState.side[state.Right].lastQRank = nextSymb;
                newState.LInfo = OnMatchL;
                newState.RInfo = OnMatchR;
                newState.NextPos = true;
                auto f = search_next_pos(newState);
                if (f) return true;
            }

            for (uint64_t i{FirstSymb}; i < Sigma; ++i) {
                auto newState = state;
                newState.e = state.e+1;
                newState.cur  = cursors[i];
                newState.side[state.Right].lastRank = i;
                if (Deletion) {
                    newState.LInfo = OnDeletionL;
                    newState.RInfo = OnDeletionR;
                    newState.NextPos = false;
                    auto f = search_next_pos(newState); // deletion occurred in query
                    if (f) return true;
                }
                if (!substitutionAllowed) continue;
                if (i == nextSymb) continue;

                newState.side[state.Right].lastQRank = nextSymb;
                newState.LInfo = OnSubstituteL;
                newState.RInfo = OnSubstituteR;
                newState.NextPos = true;
                auto f = search_next_pos(newState);
                if (f) return true;
            }

            if (Insertion) {
                if (insertionAllowed) {
                    auto newState = state;
                    newState.e = state.e+1;
                    newState.side[state.Right].lastQRank = nextSymb;
                    newState.LInfo = OnInsertionL;
                    newState.RInfo = OnInsertionR;
                    newState.NextPos = true;
                    auto f = search_next_pos(newState); // insertion occurred in query
                    if (f) return true;
                }
            }
        } else if (matchAllowed) {
            auto f = search_next_dir_no_errors(state);
            if (f) return true;
        }
        return false;
    }
    bool search_next_dir_no_errors(State state) const {
        auto loops = state.partitionEntryValue;
        auto nextSymb = decltype(query[0]){};
        for (size_t i{0}; i < loops; ++i) {
            nextSymb = query[state.Right?(state.queryPosR+i):(state.queryPosL-i)];
            state.cur = extend(state, nextSymb);
            if (state.cur.count() == 0) return false;
        }

        state.side[state.Right].lastRank = nextSymb;
        state.side[state.Right].lastQRank = nextSymb;
        state.part += 1;
        state.partitionEntryValue = 0;
        if (state.part != partition.size()) {
            state.partitionEntryValue = partition[search.pi[state.part]];
        }
        if (state.Right) {
            state.queryPosR += loops;
            state.RInfo = 'M';
        } else {
            state.queryPosL -= loops;
            state.LInfo = 'M';
        }
        bool res = search_next(state);
        return res;
    }
    bool search_next_dir_single(State const& state) const {
        char const TInfo = state.Right ? state.RInfo : state.LInfo;

        bool const Deletion     = (TInfo != 'S' && TInfo != 'I') && Edit;
        bool const Insertion    = (TInfo != 'S' && TInfo != 'D') && Edit;

        char const OnMatchL      = state.Right ? state.LInfo : 'M';
        char const OnMatchR      = state.Right ? 'M'   : state.RInfo;
        char const OnSubstituteL = state.Right ? state.LInfo : 'S';
        char const OnSubstituteR = state.Right ? 'S'   : state.RInfo;
        char const OnDeletionL   = state.Right ? state.LInfo : 'D';
        char const OnDeletionR   = state.Right ? 'D'   : state.RInfo;
        char const OnInsertionL  = state.Right ? state.LInfo : 'I';
        char const OnInsertionR  = state.Right ? 'I'   : state.RInfo;


        auto [curISymb, icursorNext] = [&]() -> std::tuple<size_t, cursor_t> {
            if (state.Right) {
                auto symb = state.cur.symbolRight();
                auto cur_  = state.cur.extendRight(symb);
                return {symb, cur_};
            } else {
                auto symb = state.cur.symbolLeft();
                auto cur_  = state.cur.extendLeft(symb);
                return {symb, cur_};
            }
        }();

        auto curQSymb = query[state.Right?state.queryPosR:state.queryPosL];

        bool insertionAllowed    = (state.partitionEntryValue > 1 or search.l[state.part] <= state.e+1)
                                   and state.e+1 <= search.u[state.part];
        bool substitutionAllowed = insertionAllowed;
        bool mismatchAllowed     = state.e+1 <= search.u[state.part];

        if (Insertion) {
            if (insertionAllowed) {
                auto newState = state;
                newState.e = state.e+1;
                newState.side[state.Right].lastQRank = curQSymb;
                newState.LInfo = OnInsertionL;
                newState.RInfo = OnInsertionR;
                newState.NextPos = true;
                bool f = search_next_pos(newState);
                if (f) return true;
            }
        }

        // only insertions are possible
        if (curISymb < FirstSymb) {
            return false;
        }

        bool matchAllowed    = (state.partitionEntryValue > 1 or search.l[state.part] <= state.e)
                               and state.e <= search.u[state.part]
                               and (TInfo != 'I' or curQSymb != state.side[state.Right].lastQRank)
                               and (TInfo != 'D' or curQSymb != state.side[state.Right].lastRank);

        if (curISymb == curQSymb) {
            if (matchAllowed) {
                if (!mismatchAllowed) {
                    auto f = search_next_dir_no_errors(state);
                    if (f) return true;
                    return false;
                }
                auto newState = state;
                newState.side[state.Right].lastRank = curQSymb;
                newState.side[state.Right].lastQRank = curQSymb;
                newState.cur = icursorNext;
                newState.LInfo = OnMatchL;
                newState.RInfo = OnMatchR;
                newState.NextPos = true;
                bool f = search_next_pos(newState);
                if (f) return true;
            }
            if (Deletion) {
                if (mismatchAllowed) {
                    auto newState = state;
                    newState.e = state.e+1;
                    newState.side[state.Right].lastRank = curISymb;
                    newState.cur = icursorNext;
                    newState.LInfo = OnDeletionL;
                    newState.RInfo = OnDeletionR;
                    newState.NextPos = false;
                    bool f = search_next_pos(newState);
                    if (f) return true;
                }
            }
        } else if (mismatchAllowed) {
            auto newState = state;
            newState.e = state.e+1;
            newState.side[state.Right].lastRank = curISymb;
            newState.cur = icursorNext;


            // search substitute
            if (substitutionAllowed) {
                auto s2 = Restore{newState.side[state.Right].lastQRank, curQSymb};
                newState.LInfo = OnSubstituteL;
                newState.RInfo = OnSubstituteR;
                newState.NextPos = true;
                bool f = search_next_pos(newState);
                if (f) return true;
            }

            if (Deletion) {
                newState.LInfo = OnDeletionL;
                newState.RInfo = OnDeletionR;
                newState.NextPos = false;
                bool f = search_next_pos(newState);
                if (f) return true;
            }
        }
        return false;
    }
};


template <bool Edit, typename index_t, Sequence query_t, typename delegate_t>
void search_impl(index_t const& index, query_t const& query, search_scheme::Scheme const& search_scheme, std::vector<size_t> const& partition, delegate_t&& delegate) {
    using cursor_t = select_cursor_t<index_t>;
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

    for (auto const& search : search_scheme) {
        bool f = Search<Edit, index_t, query_t, decltype(search), decltype(internal_delegate)>{index, query, search, partition, internal_delegate}.run();
        if (f) {
            return;
        }
    }
}

/* search until exactly n results are being found
 *
 * \tparam Edit: should Edit distance (or hamming distance be applied while searching)
 * \param index: An Index, can be FmIndex or a BiFMIndex
 * \param queries: a range of queries, do not have to be of same length
 * \param maxErrors: number of errors allowed (important to edit distance or hamming distance)
 * \param n: number of results that should be found
 * \param delegate_t: callback function to report the results, Must accept there parameters: size_t qidx, auto cur, size_t e
 *          qidx: the position of the query inside the queries range
 *          cur:  cursor of the FMIndex marking the result positions
 *          e:    number of errors that these matches have
 * \param searchSelectSearchScheme_t: callback that helps selecting a proper search scheme, Must accept one parameters: size_t length
 *          length: length of query
 */
template <bool Edit, typename index_t, Sequences queries_t, typename selectSearchScheme_t, typename delegate_t>
void search_n_impl(index_t const& index, queries_t&& queries, selectSearchScheme_t&& selectSearchScheme, delegate_t&& delegate, size_t n) {
    if (queries.empty()) return;
    if (n == 0) return;
    for (size_t qidx{}; qidx < queries.size(); ++qidx) {
        size_t ct{};
        auto const& [search_scheme, partition] = selectSearchScheme(queries[qidx].size());
        search_impl<Edit>(index, queries[qidx], search_scheme, partition, [&] (auto cur, size_t e) {
            if (cur.count() + ct > n) {
                cur.len = n-ct;
            }
            ct += cur.count();
            delegate(qidx, cur, e);
            return ct == n;
        });
    }
}


// convenience function, with passed search scheme and multiple queries
template <bool Edit=true, typename index_t, Sequences queries_t, typename delegate_t>
void search(index_t const& index, queries_t&& queries, search_scheme::Scheme const& search_scheme, std::vector<size_t> const& partition, delegate_t&& delegate, size_t n = std::numeric_limits<size_t>::max()) {
    // function that selects a search scheme
    auto selectSearchScheme = [&]([[maybe_unused]] size_t length) -> auto {
        return std::tie(search_scheme, partition);
    };
    search_n_impl<Edit>(index, queries, selectSearchScheme, delegate, n);
}

// convenience function, with auto selected search scheme and multiple queries
template <bool Edit=true, typename index_t, Sequences queries_t, typename delegate_t>
void search(index_t const& index, queries_t&& queries, size_t maxErrors, delegate_t&& delegate, size_t n = std::numeric_limits<size_t>::max()) {
    auto selectSearchScheme = [&]([[maybe_unused]] size_t length) -> auto {
        auto const& search_scheme = getCachedSearchScheme<Edit>(0, maxErrors);
        auto const& partition     = getCachedPartition(search_scheme[0].pi.size(), length);
        return std::tie(search_scheme, partition);
    };
    search_n_impl<Edit>(index, queries, selectSearchScheme, delegate, n);
}


// convenience function, with passed search scheme and multiple queries
template <bool Edit=true, typename index_t, Sequences queries_t, typename delegate_t>
void search_best(index_t const& index, queries_t&& queries, std::vector<std::tuple<search_scheme::Scheme, std::vector<size_t>>> const& search_schemes, delegate_t&& delegate, size_t n = std::numeric_limits<size_t>::max()) {
    for (size_t i{}; i < search_schemes.size(); ++i) {
        bool found = false;
        auto report = [&](size_t qidx, auto cur, size_t e) {
            found = true;
            delegate(qidx, cur, e);
        };
        auto const& [search_scheme, partition] = search_schemes[i];
        search(index, queries, search_scheme, partition, report, n);
        if (found) {
            break;
        }
    }
}

// convenience function, with auto selected search scheme and multiple queries
template <bool Edit=true, typename index_t, Sequences queries_t, typename delegate_t>
void search_best(index_t const& index, queries_t&& queries, size_t maxErrors, delegate_t&& delegate, size_t n = std::numeric_limits<size_t>::max()) {
    for (size_t i{}; i < maxErrors; ++i) {
        bool found = false;
        auto report = [&](size_t qidx, auto cur, size_t e) {
            found = true;
            delegate(qidx, cur, e);
        };
        search(index, queries, i, report, n);
        if (found) {
            break;
        }
    }
}

}
