// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "CachedSearchScheme.h"
#include "Restore.h"
#include "SelectCursor.h"

#include <array>
#include <cstddef>

/**
 * like search_ng24 but:
 *  - uses a partitioning system instead of a change of search scheme themself
 */
namespace fmc::search_ng25 {

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
    mutable std::vector<size_t> partition;
    delegate_t const& delegate;

    struct Side {
        uint8_t lastRank{};
        uint8_t lastQRank{};
    };
    mutable std::array<Side, 2> side;
    mutable size_t e{};
    mutable size_t part{};
    mutable size_t queryPosL{}, queryPosR{};

    Search(index_t const& _index, query_t const& _query, search_t const& _search, std::vector<size_t> const& _partition, delegate_t const& _delegate)
        : index     {_index}
        , query     {_query}
        , search    {_search}
        , partition {_partition}
        , delegate  {_delegate}
    {
        // check how many characters are before the first query char
        for (size_t i{0}; i < search.pi[0]; ++i) {
            queryPosL += partition[i];
            queryPosR += partition[i];
        }
        // Notice: queryPosL might 'underflow' (also at other spots)
        // but that is fine, after an underflow this variable is not
        // being accessed before overflowing again.
        queryPosL -= 1;
    }

    bool run() {
        auto cur       = cursor_t{index};
        return search_next<'M', 'M'>(cur);
    }


    template <bool Right>
    static auto extend(cursor_t const& cur, uint64_t symb) noexcept {
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
    bool search_next(cursor_t const& cur) const {
        if (cur.count() == 0) {
            return false;
        }
        if (part == partition.size()) {
            if constexpr (!Edit || ((LInfo == 'M' or LInfo == 'I') and (RInfo == 'M' or RInfo == 'I'))) {
                if (search.l.back() <= e and e <= search.u.back()) {
                    return delegate(cur, e);
                }
            }
            return false;
        }

        bool Right = (part==0) || (search.pi[part-1] < search.pi[part]);
        if (cur.count() > 1) {
            if (Right) {
                return search_next_dir<LInfo, RInfo, true>(cur);
            } else {
                return search_next_dir<LInfo, RInfo, false>(cur);
            }
        } else {
            if (Right) {
                return search_next_dir_single<LInfo, RInfo, true>(cur);
            } else {
                return search_next_dir_single<LInfo, RInfo, false>(cur);
            }
        }
    }

    template <char LInfo, char RInfo, bool Right, bool NextPos>
    bool search_next_pos(cursor_t const& cur) const {
        if (cur.count() == 0) return false;
        if constexpr (NextPos) {
            if constexpr (Right) queryPosR += 1;
            else queryPosL -= 1;
            partition[search.pi[part]] -= 1;

            if (partition[search.pi[part]] == 0) {
                part += 1;
                bool res = search_next<LInfo, RInfo>(cur);
                part -= 1;
                partition[search.pi[part]] += 1;
                if constexpr (Right) queryPosR -= 1;
                else queryPosL += 1;
                return res;
            }
        }

        bool r = [&]() {
            if (cur.count() > 1) {
                return search_next_dir<LInfo, RInfo, Right>(cur);
            } else {
                return search_next_dir_single<LInfo, RInfo, Right>(cur);
            }
        }();
        if constexpr (NextPos) {
            if constexpr (Right) queryPosR -= 1;
            else queryPosL += 1;
            partition[search.pi[part]] += 1;
        }
        return r;
    }

    template <char LInfo, char RInfo, bool Right>
    bool search_next_dir(cursor_t const& cur) const {
        static constexpr char TInfo = Right ? RInfo : LInfo;

        constexpr bool Deletion     = (TInfo != 'S' && TInfo != 'I') && Edit;
        constexpr bool Insertion    = (TInfo != 'S' && TInfo != 'D') && Edit;

        constexpr char OnMatchL      = Right ? LInfo : 'M';
        constexpr char OnMatchR      = Right ? 'M'   : RInfo;
        constexpr char OnSubstituteL = Right ? LInfo : 'S';
        constexpr char OnSubstituteR = Right ? 'S'   : RInfo;
        constexpr char OnDeletionL   = Right ? LInfo : 'D';
        constexpr char OnDeletionR   = Right ? 'D'   : RInfo;
        constexpr char OnInsertionL  = Right ? LInfo : 'I';
        constexpr char OnInsertionR  = Right ? 'I'   : RInfo;

        auto nextSymb = query[Right?queryPosR:queryPosL];

        bool matchAllowed    = (partition[search.pi[part]] > 1 or search.l[part] <= e)
                               and e <= search.u[part]
                               and (TInfo != 'I' or nextSymb != side[Right].lastQRank)
                               and (TInfo != 'D' or nextSymb != side[Right].lastRank);
        bool insertionAllowed    = (partition[search.pi[part]] > 1 or search.l[part] <= e+1)
                                   and e+1 <= search.u[part];
        bool substitutionAllowed = insertionAllowed;
        bool mismatchAllowed     = e+1 <= search.u[part];

        if (mismatchAllowed) {
            auto cursors = extend<Right>(cur);

            if (matchAllowed) {
                auto const& newCur = cursors[nextSymb];
                auto s1 = Restore{side[Right].lastRank, nextSymb};
                auto s2 = Restore{side[Right].lastQRank, nextSymb};
                auto f = search_next_pos<OnMatchL, OnMatchR, Right, /*.nextPos=*/true>(newCur);
                if (f) return true;
            }

            auto se = Restore{e, e+1};
            for (uint64_t i{FirstSymb}; i < Sigma; ++i) {
                auto const& newCur = cursors[i];

                auto s1 = Restore{side[Right].lastRank, i};
                if constexpr (Deletion) {
                    auto f = search_next_pos<OnDeletionL, OnDeletionR, Right, /*.nextPos=*/false>(newCur); // deletion occurred in query
                    if (f) return true;
                }
                if (!substitutionAllowed) continue;
                if (i == nextSymb) continue;

                auto s2 = Restore{side[Right].lastQRank, nextSymb};
                auto f = search_next_pos<OnSubstituteL, OnSubstituteR, Right, /*.nextPos=*/true>(newCur);
                if (f) return true;
            }

            if constexpr (Insertion) {
                if (insertionAllowed) {
                    auto s2 = Restore{side[Right].lastQRank, nextSymb};
                    auto f = search_next_pos<OnInsertionL, OnInsertionR, Right, /*.nextPos=*/true>(cur); // insertion occurred in query
                    if (f) return true;
                }
            }
        } else if (matchAllowed) {
            auto f = search_next_dir_no_errors<OnMatchL, OnMatchR, Right>(cur);
            if (f) return true;
        }
        return false;
    }

    template <char LInfo, char RInfo, bool Right>
    bool search_next_dir_no_errors(cursor_t cur) const {
        auto loops = partition[search.pi[part]];
        auto nextSymb = decltype(query[0]){};
        for (size_t i{0}; i < loops; ++i) {
            nextSymb = query[Right?(queryPosR+i):(queryPosL-i)];
            cur = extend<Right>(cur, nextSymb);
            if (cur.count() == 0) return false; //!TODO reset values
        }

        auto s1 = Restore{side[Right].lastRank, nextSymb};
        auto s2 = Restore{side[Right].lastQRank, nextSymb};
        partition[search.pi[part]] = 0;
        part += 1;
        if constexpr (Right) queryPosR += loops;
        else queryPosL -= loops;
        bool res = search_next<LInfo, RInfo>(cur);
        part -= 1;
        partition[search.pi[part]] = loops;
        if constexpr (Right) queryPosR -= loops;
        else queryPosL += loops;
        return res;
    }

    template <char LInfo, char RInfo, bool Right>
    bool search_next_dir_single(cursor_t const& cur) const {
        static constexpr char TInfo = Right ? RInfo : LInfo;

        constexpr bool Deletion     = (TInfo != 'S' && TInfo != 'I') && Edit;
        constexpr bool Insertion    = (TInfo != 'S' && TInfo != 'D') && Edit;

        constexpr char OnMatchL      = Right ? LInfo : 'M';
        constexpr char OnMatchR      = Right ? 'M'   : RInfo;
        constexpr char OnSubstituteL = Right ? LInfo : 'S';
        constexpr char OnSubstituteR = Right ? 'S'   : RInfo;
        constexpr char OnDeletionL   = Right ? LInfo : 'D';
        constexpr char OnDeletionR   = Right ? 'D'   : RInfo;
        constexpr char OnInsertionL  = Right ? LInfo : 'I';
        constexpr char OnInsertionR  = Right ? 'I'   : RInfo;


        auto [curISymb, icursorNext] = [&]() -> std::tuple<size_t, cursor_t> {
            if constexpr (Right) {
                auto symb = cur.symbolRight();
                auto cur_  = cur.extendRight(symb);
                return {symb, cur_};
            } else {
                auto symb = cur.symbolLeft();
                auto cur_  = cur.extendLeft(symb);
                return {symb, cur_};
            }
        }();

        auto curQSymb = query[Right?queryPosR:queryPosL];

        bool insertionAllowed    = search.l[part] <= e+partition[search.pi[part]] and e+1 <= search.u[part];
        bool substitutionAllowed = insertionAllowed;
        bool mismatchAllowed     = e+1 <= search.u[part];

        if constexpr (Insertion) {
            if (insertionAllowed) {
                auto se = Restore{e, e+1};
                auto s2 = Restore{side[Right].lastQRank, curQSymb};

                bool f = search_next_pos<OnInsertionL, OnInsertionR, Right, /*.nextPos=*/true>(cur);
                if (f) return true;
            }
        }

        // only insertions are possible
        if (curISymb < FirstSymb) {
            return false;
        }

        bool matchAllowed    = (partition[search.pi[part]] > 1 or search.l[part] <= e)
                               and e <= search.u[part]
                               and (TInfo != 'I' or curQSymb != side[Right].lastQRank)
                               and (TInfo != 'D' or curQSymb != side[Right].lastRank);

        if (curISymb == curQSymb) {
            if (matchAllowed) {
                if (!mismatchAllowed) {
                    auto f = search_next_dir_no_errors<OnMatchL, OnMatchR, Right>(cur);
                    if (f) return true;
                    return false;
                }
                auto s1 = Restore{side[Right].lastRank, curQSymb};
                auto s2 = Restore{side[Right].lastQRank, curQSymb};
                bool f = search_next_pos<OnMatchL, OnMatchR, Right, /*.nextPos=*/true>(icursorNext);
                if (f) return true;
            }
            if constexpr (Deletion) {
                if (mismatchAllowed) {
                    auto se = Restore{e, e+1};
                    auto s1 = Restore{side[Right].lastRank, curISymb};
                    bool f = search_next_pos<OnDeletionL, OnDeletionR, Right, /*.nextPos=*/false>(icursorNext);
                    if (f) return true;
                }
            }
        } else if (mismatchAllowed) {
            auto se = Restore{e, e+1};

            // search substitute
            if (substitutionAllowed) {
                auto s1 = Restore{side[Right].lastRank, curISymb};
                auto s2 = Restore{side[Right].lastQRank, curQSymb};

                bool f = search_next_pos<OnSubstituteL, OnSubstituteR, Right, /*.nextPos=*/true>(icursorNext);
                if (f) return true;
            }

            if constexpr (Deletion) {
                auto s1 = Restore{side[Right].lastRank, curISymb};
                bool f = search_next_pos<OnDeletionL, OnDeletionR, Right, /*.nextPos=*/false>(icursorNext);
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
void search_n_impl(index_t const& index, queries_t&& queries, selectSearchScheme_t&& selectSearchScheme, std::vector<size_t> const& partition, delegate_t&& delegate, size_t n) {
    if (queries.empty()) return;
    if (n == 0) return;
    for (size_t qidx{}; qidx < queries.size(); ++qidx) {
        size_t ct{};
        auto const& search_scheme = selectSearchScheme(queries[qidx].size());
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
    auto selectSearchScheme = [&]([[maybe_unused]] size_t length) -> auto& {
        return search_scheme;
    };
    search_n_impl<Edit>(index, queries, selectSearchScheme, partition, delegate, n);
}

}
