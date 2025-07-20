// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "SelectCursor.h"

#include <array>
#include <cstddef>

/**
 * like search_ng21 but:
 *  - correct handling of left/right avoidance of ID and other combinations
 *  - handling of matches followed by insertions or deletions
 */
namespace fmc::search_ng24 {

template <typename T, typename T2>
struct Store {
    T* value;
    T oldValue;
    Store(T& _value, T2 newValue)
        : value{&_value}
        , oldValue{*value}
    {
        *value = static_cast<T>(newValue);
    }
    Store(Store&& store) = delete;
    ~Store() {
        *value = oldValue;
    }
};


template <typename index_t, typename query_t, typename search_t, typename delegate_t>
struct Search {
    constexpr static size_t Sigma = index_t::Sigma;

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
            if constexpr ((LInfo == 'M' or LInfo == 'I') and (RInfo == 'M' or RInfo == 'I')) {
                return delegate(cur, e);
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

        constexpr bool Deletion     = TInfo != 'S' && TInfo != 'I';
        constexpr bool Insertion    = TInfo != 'S' && TInfo != 'D';

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
                auto s1 = Store{side[Right].lastRank, nextSymb};
                auto s2 = Store{side[Right].lastQRank, nextSymb};
                auto sp = Store{part, part+1};
                bool f = search_next<OnMatchL, OnMatchR>(newCur);
                if (f) return true;
            }

            auto se = Store{e, e+1};
            for (uint64_t i{1}; i < Sigma; ++i) {
                auto const& newCur = cursors[i];

                auto s1 = Store{side[Right].lastRank, i};
                if constexpr (Deletion) {
                    if (TInfo != 'M' || side[Right].lastQRank != i) {
                        bool f = search_next<OnDeletionL, OnDeletionR>(newCur); // deletion occurred in query
                        if (f) return true;
                    }
                }

                if (i == nextSymb) continue;

                auto s2 = Store{side[Right].lastQRank, nextSymb};
                auto sp = Store{part, part+1};
                bool f = search_next<OnSubstituteL, OnSubstituteR>(newCur); // as substitution
                if (f) return true;
            }

            if constexpr (Insertion) {
                if (TInfo != 'M' || side[Right].lastQRank != nextSymb) {
                    auto s2 = Store{side[Right].lastQRank, nextSymb};
                    auto sp = Store{part, part+1};
                    bool f = search_next<OnInsertionL, OnInsertionR>(cur); // insertion occurred in query
                    if (f) return true;
                }
            }
        } else if (matchAllowed) {
            auto s1 = Store{side[Right].lastRank, nextSymb};
            auto s2 = Store{side[Right].lastQRank, nextSymb};
            auto newCur = extend<Right>(cur, nextSymb);
            auto sp = Store{part, part+1};
            bool f = search_next<OnMatchL, OnMatchR>(newCur);
            if (f) return true;
        }
        return false;
    }

    template <char LInfo, char RInfo, bool Right>
    bool search_next_dir_single(cursor_t const& cur) const {
        static constexpr char TInfo = Right ? RInfo : LInfo;

        constexpr bool Deletion     = TInfo != 'S' && TInfo != 'I';
        constexpr bool Insertion    = TInfo != 'S' && TInfo != 'D';

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
                auto symb = index.bwtRev.symbol(cur.lbRev);
                auto cur_  = cur.extendRight(symb);
                return {symb, cur_};
            } else {
                auto symb = index.bwt.symbol(cur.lb);
                auto cur_  = cur.extendLeft(symb);
                return {symb, cur_};
            }
        }();

        auto curQSymb = query[search.pi[part]];

        bool mismatchAllowed = search.l[part] <= e+1 and e+1 <= search.u[part];


        if constexpr (Insertion) {
            if (mismatchAllowed) {
                if (TInfo != 'M' || side[Right].lastQRank != curQSymb) {
                    auto se = Store{e, e+1};
                    auto s2 = Store{side[Right].lastQRank, curQSymb};
                    auto sp = Store{part, part+1};

                    bool f = search_next<OnInsertionL, OnInsertionR>(cur);
                    if (f) return true;
                }
            }
        }




        // only insertions are possible
        if (curISymb == 0) {
            return false;
        }


        bool matchAllowed    = search.l[part] <= e and e <= search.u[part]
                               and (TInfo != 'I' or curQSymb != side[Right].lastQRank)
                               and (TInfo != 'D' or curQSymb != side[Right].lastRank);

        if (curISymb == curQSymb) {
            if (matchAllowed) {
                auto s1 = Store{side[Right].lastRank, curQSymb};
                auto s2 = Store{side[Right].lastQRank, curQSymb};
                auto sp = Store{part, part+1};
                bool f = search_next<OnMatchL, OnMatchR>(icursorNext);
                if (f) return true;
            }
            if constexpr (Deletion) {
                if (mismatchAllowed) {
                    if (TInfo != 'M' || side[Right].lastQRank != curISymb) {
                        auto se = Store{e, e+1};
                        auto s1 = Store{side[Right].lastRank, curISymb};
                        bool f = search_next<OnDeletionL, OnDeletionR>(icursorNext);
                        if (f) return true;
                    }
                }
            }
        } else if (mismatchAllowed) {
            auto se = Store{e, e+1};

            // search substitute
            {
                auto s1 = Store{side[Right].lastRank, curISymb};
                auto s2 = Store{side[Right].lastQRank, curQSymb};
                auto sp = Store{part, part+1};

                bool f = search_next<OnSubstituteL, OnSubstituteR>(icursorNext);
                if (f) return true;
            }

            if constexpr (Deletion) {
                if (TInfo != 'M' || side[Right].lastQRank != curISymb) {
                    auto s1 = Store{side[Right].lastRank, curISymb};
                    bool f = search_next<OnDeletionL, OnDeletionR>(icursorNext);
                    if (f) return true;
                }
            }
        }
        return false;
    }

};



template <typename index_t, typename query_t, typename search_scheme_t, typename delegate_t>
void search_with_abort(index_t const& index, query_t const& query, search_scheme_t const& search_scheme, delegate_t&& delegate) {
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
        bool f = Search{index, query, search, internal_delegate}.run();
        if (f) {
            return;
        }
    }
}

template <typename index_t, typename queries_t, typename search_scheme_t, typename delegate_t>
void search(index_t const& index, queries_t&& queries, search_scheme_t const& search_scheme, delegate_t&& delegate) {
    if (search_scheme.empty()) return;

    for (size_t qidx{}; qidx < queries.size(); ++qidx) {
        search_with_abort(index, queries[qidx], search_scheme, [&](auto const& cur, size_t e) {
            delegate(qidx, cur, e);
        });
    }
}

template <typename index_t, typename queries_t, typename search_scheme_t, typename delegate_t>
void search_n(index_t const& index, queries_t&& queries, search_scheme_t const& search_scheme, size_t n, delegate_t&& delegate) {
    if (search_scheme.empty()) return;

    for (size_t qidx{}; qidx < queries.size(); ++qidx) {
        size_t ct{};
        search_with_abort(index, queries[qidx], search_scheme, [&] (auto cur, size_t e) {
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
