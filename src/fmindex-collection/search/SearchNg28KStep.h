// SPDX-FileCopyrightText: 2026 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "CachedSearchScheme.h"
#include "Restore.h"
#include "SelectCursor.h"

#include <array>
#include <cstddef>

/**
 * like search_ng25 but:
 *  - optimized further reduced templates
 *  - optimized 'Insert' case
 */
namespace fmc::search_ng28_kstep {

template <bool Edit, typename index_t, typename query_t, typename search_t, typename delegate_t>
struct Search {
    using cursor_t = select_cursor_t<index_t>;

    constexpr static size_t Sigma = index_t::Sigma;
    constexpr static size_t FirstSymb = []() -> size_t {
        if constexpr (requires() { { index_t::FirstSymb }; }) {
            return index_t::FirstSymb;
        }
        return 1;
    }();

    constexpr bool static HasKStep = requires() {
        { index_t::KStep };
    };
    constexpr size_t static KStep = []() constexpr -> size_t { if constexpr(HasKStep) return index_t::KStep; return 1ul; }();
    constexpr size_t static SigmaKStep = []() constexpr -> size_t { if constexpr(HasKStep) return index_t::SigmaKStep; return Sigma; }();


    index_t const& index;
    query_t const& query;
    search_t const& search;
    std::vector<size_t> partition;
    delegate_t const& delegate;

    enum dir_t : int {
        Left = -1,
        Right = 1,
    };

    struct Side {
        uint8_t lastRank{};
        uint8_t lastQRank{};
        char    info{'M'};
        size_t  queryPos{};
    };
    mutable std::array<Side, 3> sideImpl;
    mutable Side* side{nullptr};
    mutable size_t e{};
    mutable size_t part{std::numeric_limits<size_t>::max()}; // first action will be overflow, this is desired effect
    mutable dir_t dir{dir_t::Right};
    mutable size_t partitionPart{};
    mutable uint8_t lb, ub;

    Search(index_t const& _index, query_t const& _query, search_t const& _search, std::vector<size_t> const& _partition, delegate_t const& _delegate)
        : index     {_index}
        , query     {_query}
        , search    {_search}
        , partition {_partition}
        , delegate  {_delegate}
    {
        // check how many characters are before the first query char
        for (size_t i{0}; i < search.pi[0]; ++i) {
            sideImpl[0].queryPos += partition[i];
            sideImpl[2].queryPos += partition[i];
        }
        // Notice: sideImpl[0].queryPos might 'underflow' (also at other spots)
        // but that is fine, after an underflow this variable is not
        // being accessed before overflowing again.
        sideImpl[0].queryPos -= 1;
    }

    bool run() {
        auto cur = cursor_t{index};
        return search_next_part(cur);
    }

    auto extend(cursor_t const& cur, uint64_t symb) const noexcept {
        if (dir == dir_t::Right) return cur.extendRight(symb);
        return cur.extendLeft(symb);
    }

    auto extend(cursor_t const& cur) const noexcept {
        if (dir == dir_t::Right) return cur.extendRight();
        return cur.extendLeft();
    }

    auto computeErrorCases() const {
        auto const tinfo               = side->info;
        bool const mismatchAllowed     = e+1 <= ub;
        if (!mismatchAllowed) return std::make_tuple(false, false, false, false);
        bool const substitutionAllowed = ((Edit && (partitionPart > 1 || lb <= e+1)) || (!Edit && lb <= e+partitionPart)) && mismatchAllowed;
        bool const insertionAllowed    = Edit && tinfo != 'S' && tinfo != 'D' && substitutionAllowed;
        bool const deletionAllowed     = Edit && tinfo != 'S' && tinfo != 'I' && mismatchAllowed;
        return std::make_tuple(mismatchAllowed, substitutionAllowed, insertionAllowed, deletionAllowed);
    }
    auto computeMatchCase(size_t nextSymb) const {
        auto const tinfo        = side->info;
        bool const matchAllowed = ((!Edit && (partitionPart > 1 || lb <= e))
                                    || (Edit && lb <= e + partitionPart-1))
                                 and e <= ub
                                 and (tinfo != 'I' or nextSymb != side->lastQRank)
                                 and (tinfo != 'D' or nextSymb != side->lastRank);
        return matchAllowed;
    }

    /*
     * Searches for the next part
     */
    bool search_next_part(cursor_t const& cur) const {
        if (cur.count() == 0) {
            return false;
        }
        auto r_p = RestoreAdd{part, 1};
        if (part == partition.size()) {
            auto il = sideImpl[0].info;
            auto ir = sideImpl[2].info;
            if (!Edit || ((il == 'M' or il == 'I') and (ir == 'M' or ir == 'I'))) {
                if (search.l.back() <= e and e <= search.u.back()) {
                    return delegate(cur, e);
                }
            }
            return false;
        }

        bool right = (part==0) || (search.pi[part-1] < search.pi[part]);
        auto r_dir = Restore{dir, right?dir_t::Right:dir_t::Left};
        auto r_side = Restore{side, right?&sideImpl[2]:&sideImpl[0]};
        auto r_part = Restore{partitionPart, partition[search.pi[part]]};
        auto r_lb   = Restore{lb, search.l[part]};
        auto r_ub   = Restore{ub, search.u[part]};

        return search_next_symb(cur);
    }

    /**
     * move search position to next character, check if partition is ending
     * if required start a new part/partition search
     */
    bool check_and_search_next_symb(cursor_t const& cur) const {
        if (cur.count() == 0) return false;
        auto r_qp = RestoreAdd{side->queryPos, dir};
        auto r_p  = RestoreSub{partitionPart, 1};

        if (partitionPart == 0) {
            return search_next_part(cur);
        }
        return search_next_symb(cur);
    }

    /* match next symbol
     */
    bool search_next_symb(cursor_t const& cur) const {
        if (cur.count() == 1) {
            return search_next_symb_single(cur);
        }

        auto nextSymb = query[side->queryPos];

        auto matchAllowed = computeMatchCase(nextSymb);
        auto [mismatchAllowed, substitutionAllowed, insertionAllowed, deletionAllowed] = computeErrorCases();

        if (!mismatchAllowed && matchAllowed) {
            return search_next_part_no_errors(cur);
        }

        if (mismatchAllowed) {
            auto r_i   = Restore{side->info};
            auto r_lqr = Restore{side->lastQRank};
            auto r_qp  = Restore{side->queryPos};
            auto r_p   = Restore{partitionPart};
            auto r_e   = Restore{e};

            auto cursors = extend(cur);

            while (true) {
                if (matchAllowed) {
                    auto const& newCur = cursors[nextSymb];
                    auto r_lr  = Restore{side->lastRank, nextSymb};
                    auto r_lqr = Restore{side->lastQRank, nextSymb};
                    side->info = 'M';
                    auto f = check_and_search_next_symb(newCur);
                    if (f) return true;
                }
                if (!mismatchAllowed) break;
                e += 1;
                for (uint64_t i{FirstSymb}; i < Sigma; ++i) {
                    auto const& newCur = cursors[i];
                    if (newCur.count() == 0) continue;

                    auto r_lr = Restore{side->lastRank, i};
                    if (deletionAllowed) {
                        side->info = 'D';
                        auto f = search_next_symb(newCur); // deletion occurred in query
                        if (f) return true;
                    }
                    if (!substitutionAllowed) continue;
                    if (i == nextSymb) continue;

                    auto r_lqr = Restore{side->lastQRank, nextSymb};
                    side->info = 'S';
                    auto f = check_and_search_next_symb(newCur);
                    if (f) return true;
                }

                if (!insertionAllowed) break;
                side->lastQRank = nextSymb;
                side->info = 'I';
                side->queryPos += dir;
                partitionPart -= 1;
                if (partitionPart == 0) {
                    return search_next_part(cur);
                }

                nextSymb = query[side->queryPos];

                matchAllowed = computeMatchCase(nextSymb);
                std::tie(mismatchAllowed, substitutionAllowed, insertionAllowed, deletionAllowed) = computeErrorCases();
            }
        }
        return false;
    }

    /** Searches the (rest of) a part without any errors
     *
     * May only be called if no errors are allowed
     */
    bool search_next_part_no_errors(cursor_t cur) const {
        assert(e == ub);


        auto loops = partitionPart;
        size_t i{};
        bool moveRight = (dir == dir_t::Right);
        if constexpr (HasKStep) {
            for (; i+KStep < loops; i += KStep) {
                size_t kSymb{};
                for (size_t j{}; j < KStep; ++j) {
                    auto pos = side->queryPos + (i+j)*dir ;
                    kSymb = kSymb*Sigma + query[pos];
                }
                if (cur.count() == 1 && false) {
                    auto s = moveRight?cur.symbolRightKStep():cur.symbolLeftKStep();
                    if (s != kSymb) return false;
/*                    cur = [&]() {
                        if (dir == dir_t::Right) return cur.extendRightKStep(kSymb);
                        return cur.extendLeftKStep(kSymb);
                    }();*/

                    cur = moveRight?cur.extendRightBySymbolKStep(kSymb):cur.extendLeftBySymbolKStep(kSymb);
                } else {
                    cur = moveRight?cur.extendRightKStep(kSymb):cur.extendLeftKStep(kSymb);
                }
                if (cur.count() == 0) return false; // early abort, if results are empty
            }
        }
        for (; i < loops; ++i) {
            auto pos = side->queryPos + i*dir;
            auto nextSymb = query[pos];
            if (cur.count() == 1) {
                auto s = (dir == dir_t::Right)?cur.symbolRight():cur.symbolLeft();
                if (s != nextSymb) return false;
                if constexpr (requires() { { cur.extendRightBySymbol(nextSymb) }; }) {
                    cur = (dir == dir_t::Right)?cur.extendRightBySymbol(nextSymb):cur.extendLeftBySymbol(nextSymb);
                } else {
                    cur = (dir == dir_t::Right)?cur.extendRight(nextSymb):cur.extendLeft(nextSymb);
                }
            } else {
                cur = extend(cur, nextSymb);
            }
            if (cur.count() == 0) return false; // early abort, if results are empty
        }

        auto nextSymb = query[side->queryPos + (partitionPart-1)*dir];

        auto r_lr   = Restore{side->lastRank, nextSymb};
        auto r_lqr  = Restore{side->lastQRank, nextSymb};
        auto r_p    = Restore{partitionPart, 0};
        auto r_qp   = RestoreAdd{side->queryPos, loops * dir};
        auto r_i    = Restore{side->info, 'M'};
        return search_next_part(cur);
    }

    /** Assumes only a single row is marked.
     *
     * Searches for the next symbol
     */
    bool search_next_symb_single(cursor_t const& cur) const {
        auto r_e    = Restore{e};
        auto r_lqr  = Restore{side->lastQRank};
        auto r_i    = Restore{side->info};
        auto r_qp   = Restore{side->queryPos};
        auto r_p    = Restore{partitionPart};

        // detect next symbol/cursor
        auto [curISymb, icursorNext] = [&]() -> std::tuple<size_t, cursor_t> {
            if constexpr (requires() { { cur.extendRightBySymbol() }; }) {
                if (dir == dir_t::Right) {
                    return cur.extendRightBySymbol();
                } else {
                    return cur.extendLeftBySymbol();
                }
            } else {
                if (dir == dir_t::Right) {
                    auto symb = cur.symbolRight();
                    auto cur_ = cur.extendRight(symb);
                    return {symb, cur_};
                } else {
                    auto symb = cur.symbolLeft();
                    auto cur_ = cur.extendLeft(symb);
                    return {symb, cur_};
                }
            }
        }();

        while (true) {
            auto [mismatchAllowed, substitutionAllowed, insertionAllowed, deletionAllowed] = computeErrorCases();
            auto curQSymb = query[side->queryPos];

            // only insertions are possible
            if (curISymb >= FirstSymb) {
                auto const matchAllowed = computeMatchCase(curQSymb);

                if (curISymb == curQSymb) {
                    if (matchAllowed) {
                        if (!mismatchAllowed) {
                            return search_next_part_no_errors(cur);
                        }
                        auto r_lr  = Restore{side->lastRank, curISymb};
                        auto r_lqr = Restore{side->lastQRank, curQSymb};
                        side->info = 'M';
                        bool f = check_and_search_next_symb(icursorNext);
                        if (f) return true;
                    }
                } else if (substitutionAllowed) {
                    // search substitute
                    auto r_e = Restore{e, e+1};
                    auto r_lr  = Restore{side->lastRank, curISymb};
                    auto r_lqr = Restore{side->lastQRank, curQSymb};
                    side->info = 'S';
                    bool f = check_and_search_next_symb(icursorNext);
                    if (f) return true;
                }
                if (deletionAllowed) {
                    auto r_e  = Restore{e, e+1};
                    auto r_lr = Restore{side->lastRank, curISymb};
                    side->info = 'D';
                    bool f = search_next_symb_single(icursorNext);
                    if (f) return true;
                }
            }

            if (!insertionAllowed) break;
            e += 1;
            side->lastQRank = curQSymb;
            side->info = 'I';

            side->queryPos += dir;
            partitionPart -= 1;
            if (partitionPart == 0) {
                return search_next_part(cur);
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
