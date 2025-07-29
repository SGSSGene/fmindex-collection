// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "CachedSearchScheme.h"
#include "Restore.h"
#include "SelectCursor.h"

#include <array>
#include <cstddef>

/**
 * like search_ng21 but:
 *  - correct handling of left/right avoidance of ID and other combinations
 *  - handling of matches followed by insertions or deletions
 */
namespace fmc::search_ng24 {

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
    delegate_t const& delegate;

    struct Side {
        uint8_t lastRank{};
        uint8_t lastQRank{};
    };
    mutable std::array<Side, 2> side;
    mutable size_t e{};
    mutable size_t part{};

    Search(index_t const& _index, query_t const& _query, search_t const& _search, delegate_t const& _delegate)
        : index     {_index}
        , query     {_query}
        , search    {_search}
        , delegate  {_delegate}
    {
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

        if (part == search.pi.size()) {
            if constexpr (!Edit || ((LInfo == 'M' or LInfo == 'I') and (RInfo == 'M' or RInfo == 'I'))) {
                if (search.l[part-1] <= e and e <= search.u[part-1]) {
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

        auto nextSymb = query[search.pi[part]];

        bool matchAllowed    = search.l[part] <= e and e <= search.u[part]
                               and (TInfo != 'I' or nextSymb != side[Right].lastQRank)
                               and (TInfo != 'D' or nextSymb != side[Right].lastRank);
        bool mismatchAllowed = search.l[part] <= e+1 and e+1 <= search.u[part];

        if (mismatchAllowed) {
            auto cursors = extend<Right>(cur);

            if (matchAllowed) {
                auto const& newCur = cursors[nextSymb];
                auto s1 = Restore{side[Right].lastRank, nextSymb};
                auto s2 = Restore{side[Right].lastQRank, nextSymb};
                auto sp = Restore{part, part+1};
                bool f = search_next<OnMatchL, OnMatchR>(newCur);
                if (f) return true;
            }

            auto se = Restore{e, e+1};
            for (uint64_t i{FirstSymb}; i < Sigma; ++i) {
                auto const& newCur = cursors[i];

                auto s1 = Restore{side[Right].lastRank, i};
                if constexpr (Deletion) {
                    if (TInfo != 'M' || side[Right].lastQRank != i) {
                        bool f = search_next<OnDeletionL, OnDeletionR>(newCur); // deletion occurred in query
                        if (f) return true;
                    }
                }

                if (i == nextSymb) continue;

                auto s2 = Restore{side[Right].lastQRank, nextSymb};
                auto sp = Restore{part, part+1};
                bool f = search_next<OnSubstituteL, OnSubstituteR>(newCur); // as substitution
                if (f) return true;
            }

            if constexpr (Insertion) {
                if (TInfo != 'M' || side[Right].lastQRank != nextSymb) {
                    auto s2 = Restore{side[Right].lastQRank, nextSymb};
                    auto sp = Restore{part, part+1};
                    bool f = search_next<OnInsertionL, OnInsertionR>(cur); // insertion occurred in query
                    if (f) return true;
                }
            }
        } else if (matchAllowed) {
            auto s1 = Restore{side[Right].lastRank, nextSymb};
            auto s2 = Restore{side[Right].lastQRank, nextSymb};
            auto newCur = extend<Right>(cur, nextSymb);
            auto sp = Restore{part, part+1};
            bool f = search_next<OnMatchL, OnMatchR>(newCur);
            if (f) return true;
        }
        return false;
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

        auto curQSymb = query[search.pi[part]];

        bool mismatchAllowed = search.l[part] <= e+1 and e+1 <= search.u[part];


        if constexpr (Insertion) {
            if (mismatchAllowed) {
                if (TInfo != 'M' || side[Right].lastQRank != curQSymb) {
                    auto se = Restore{e, e+1};
                    auto s2 = Restore{side[Right].lastQRank, curQSymb};
                    auto sp = Restore{part, part+1};

                    bool f = search_next<OnInsertionL, OnInsertionR>(cur);
                    if (f) return true;
                }
            }
        }

        // only insertions are possible
        if (curISymb < FirstSymb) {
            return false;
        }

        bool matchAllowed    = search.l[part] <= e and e <= search.u[part]
                               and (TInfo != 'I' or curQSymb != side[Right].lastQRank)
                               and (TInfo != 'D' or curQSymb != side[Right].lastRank);

        if (curISymb == curQSymb) {
            if (matchAllowed) {
                auto s1 = Restore{side[Right].lastRank, curQSymb};
                auto s2 = Restore{side[Right].lastQRank, curQSymb};
                auto sp = Restore{part, part+1};
                bool f = search_next<OnMatchL, OnMatchR>(icursorNext);
                if (f) return true;
            }
            if constexpr (Deletion) {
                if (mismatchAllowed) {
                    if (TInfo != 'M' || side[Right].lastQRank != curISymb) {
                        auto se = Restore{e, e+1};
                        auto s1 = Restore{side[Right].lastRank, curISymb};
                        bool f = search_next<OnDeletionL, OnDeletionR>(icursorNext);
                        if (f) return true;
                    }
                }
            }
        } else if (mismatchAllowed) {
            auto se = Restore{e, e+1};

            // search substitute
            {
                auto s1 = Restore{side[Right].lastRank, curISymb};
                auto s2 = Restore{side[Right].lastQRank, curQSymb};
                auto sp = Restore{part, part+1};

                bool f = search_next<OnSubstituteL, OnSubstituteR>(icursorNext);
                if (f) return true;
            }

            if constexpr (Deletion) {
                if (TInfo != 'M' || side[Right].lastQRank != curISymb) {
                    auto s1 = Restore{side[Right].lastRank, curISymb};
                    bool f = search_next<OnDeletionL, OnDeletionR>(icursorNext);
                    if (f) return true;
                }
            }
        }
        return false;
    }

};



template <bool Edit=true, typename index_t, Sequence query_t, typename delegate_t>
void search(index_t const& index, query_t const& query, search_scheme::Scheme const& search_scheme, delegate_t&& delegate) {
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
        bool f = Search<Edit, index_t, query_t, decltype(search), decltype(internal_delegate)>{index, query, search, internal_delegate}.run();
        if (f) {
            return;
        }
    }
}

template <bool Edit=true, typename index_t, Sequences queries_t, typename delegate_t>
void search(index_t const& index, queries_t&& queries, search_scheme::Scheme const& search_scheme, delegate_t&& delegate) {
    if (search_scheme.empty()) return;

    for (size_t qidx{}; qidx < queries.size(); ++qidx) {
        search<Edit>(index, queries[qidx], search_scheme, [&](auto const& cur, size_t e) {
            delegate(qidx, cur, e);
        });
    }
}

template <bool Edit=true, typename index_t, Sequence query_t, typename delegate_t>
void search(index_t const& index, query_t const& query, size_t maxErrors, delegate_t&& delegate) {
    auto const& search_scheme = getCachedSearchScheme<Edit>(query.size(), 0, maxErrors);
    search<Edit>(index, query, search_scheme, delegate);
}

template <bool Edit=true, typename index_t, Sequences queries_t, typename delegate_t>
void search(index_t const& index, queries_t&& queries, size_t maxErrors, delegate_t&& delegate) {
    for (size_t qidx{}; qidx < queries.size(); ++qidx) {
        auto const& search_scheme = getCachedSearchScheme<Edit>(queries[qidx].size(), 0, maxErrors);
        search<Edit>(index, queries[qidx], search_scheme, [&](auto const& cur, size_t e) {
            delegate(qidx, cur, e);
        });
    }
}


template <bool Edit=true, typename index_t, Sequence query_t, typename delegate_t>
void search_n(index_t const& index, query_t&& query, search_scheme::Scheme const& search_scheme, size_t n, delegate_t&& delegate) {
    if (search_scheme.empty()) return;

    size_t ct{};
    search<Edit>(index, query, search_scheme, [&] (auto cur, size_t e) {
        if (cur.count() + ct > n) {
            cur.len = n-ct;
        }
        ct += cur.count();
        delegate(cur, e);
        return ct == n;
    });
}

template <bool Edit=true, typename index_t, Sequences queries_t, typename delegate_t>
void search_n(index_t const& index, queries_t&& queries, search_scheme::Scheme const& search_scheme, size_t n, delegate_t&& delegate) {
    if (search_scheme.empty()) return;

    for (size_t qidx{}; qidx < queries.size(); ++qidx) {
        size_t ct{};
        search<Edit>(index, queries[qidx], search_scheme, [&] (auto cur, size_t e) {
            if (cur.count() + ct > n) {
                cur.len = n-ct;
            }
            ct += cur.count();
            delegate(qidx, cur, e);
            return ct == n;
        });
    }
}

template <bool Edit=true, typename index_t, Sequence query_t, typename delegate_t>
void search_n(index_t const& index, query_t const& query, size_t maxErrors, size_t n, delegate_t&& delegate) {
    auto const& search_scheme = getCachedSearchScheme<Edit>(query.size(), 0, maxErrors);
    search_n<Edit>(index, query, search_scheme, n, delegate);
}


template <bool Edit=true, typename index_t, Sequences queries_t, typename delegate_t>
void search_n(index_t const& index, queries_t&& queries, size_t maxErrors, size_t n, delegate_t&& delegate) {
    for (size_t qidx{}; qidx < queries.size(); ++qidx) {
        size_t ct{};
        auto const& search_scheme = getCachedSearchScheme<Edit>(queries[qidx].size(), 0, maxErrors);
        search<Edit>(index, queries[qidx], search_scheme, [&] (auto cur, size_t e) {
            if (cur.count() + ct > n) {
                cur.len = n-ct;
            }
            ct += cur.count();
            delegate(qidx, cur, e);
            return ct == n;
        });
    }
}
}
