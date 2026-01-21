// SPDX-FileCopyrightText: 2026 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "CachedSearchScheme.h"
#include "Restore.h"
#include "SelectCursor.h"

#include <array>
#include <cstddef>

/**
 * like search_ng28 but:
 *  - optimized for KStep approach
 */
namespace fmc::search_ng28_kstep {

template <bool Edit, typename index_t, typename query_t, typename search_t, typename delegate_t>
struct Search {
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
    using cursor_t = select_cursor_t<index_t>;

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
        std::vector<size_t> symbBuffer;
        size_t usedSymbBuffer{};
    };
    mutable std::array<Side, 3> sideImpl;
    mutable Side* side{nullptr};
    mutable size_t e{};
    mutable size_t part{std::numeric_limits<size_t>::max()}; // first action will be overflow, this is desired effect
    mutable dir_t dir{dir_t::Right};
    mutable size_t partitionPart{};
    mutable uint8_t lb, ub;

    mutable std::array<cursor_t, SigmaKStep> const* nextCursors{};
    mutable cursor_t const* nextSingleCursor{};
    mutable std::array<size_t, KStep> const* nextSingleCursorSymb{};

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
        auto largestAllowedError = size_t{};
        for (auto u : search.u) {
            largestAllowedError = std::max(largestAllowedError, u);
        }
        sideImpl[0].symbBuffer.resize(query.size() + largestAllowedError);
        sideImpl[2].symbBuffer.resize(query.size() + largestAllowedError);

        // Notice: sideImpl[0].queryPos might 'underflow' (also at other spots)
        // but that is fine, after an underflow this variable is not
        // being accessed before overflowing again.
        sideImpl[0].queryPos -= 1;
    }

    bool run() {
        auto cur = cursor_t{index};
        return search_next_part(cur);
    }

/*    auto flushCursor(cursor_t cur) const noexcept {
        if constexpr (HasKStep) {
            assert(sideImpl[0].usedSymbBuffer == 0 || sideImpl[2].usedSymbBuffer == 0);
            auto flush = [&](size_t side) {
                auto rest = sideImpl[side].usedSymbBuffer % KStep;
                if (rest % KStep == 0) return;
                auto buffer = std::span{sideImpl[side].symbBuffer}.first(sideImpl[side].usedSymbBuffer);
                for (auto c : buffer.last(rest)) {
                    if (side == 0) cur = cur.extendLeft(c);
                    else cur = cur.extendRight(c);
                }
            };
            flush(0);
            flush(2);
        }
        return cur;
    }*/

    auto clearCursor(cursor_t const& _cur, auto&& cb) const noexcept {
        if constexpr (HasKStep) {
            assert((sideImpl[0].usedSymbBuffer % KStep == 0) || (sideImpl[2].usedSymbBuffer % KStep == 0));
            if (auto& side = sideImpl[0]; side.usedSymbBuffer % KStep > 0) {
                auto cur = _cur;
                auto offset = side.usedSymbBuffer % KStep;
                auto start  = side.usedSymbBuffer - offset;
                for (size_t i{0}; i < offset; ++i) {
                    cur = extend(cur, side.symbBuffer[start+i]);
                }
                auto r_used = RestoreAdd{side.usedSymbBuffer, KStep - offset};
                auto f = cb(cur);
                return f;
            }

            if (auto& side = sideImpl[2]; side.usedSymbBuffer % KStep > 0) {
                auto cur = _cur;
                auto offset = side.usedSymbBuffer % KStep;
                auto start  = side.usedSymbBuffer - offset;
                for (size_t i{0}; i < offset; ++i) {
                    cur = extend(cur, side.symbBuffer[start+i]);
                }
                auto r_used = RestoreAdd{side.usedSymbBuffer, KStep - offset};
                auto f = cb(cur);
                return f;
            }
        }
        return cb(_cur);
    }


    auto extend(cursor_t const& cur, uint64_t symb) const noexcept {
        if (dir == dir_t::Right) {
            return cur.extendRight(symb);
        } else {
            return cur.extendLeft(symb);
        }
    }

    auto extendKStep(cursor_t const& cur, std::span<size_t const, KStep> symb) const noexcept {
        if constexpr (HasKStep) {
            if (dir == dir_t::Right) return cur.extendRightKStep(symb);
            else return cur.extendLeftKStep(symb);
        } else {
            static_assert(KStep == 1);
            return extend(cur, symb[0]);
        }
    }

    auto extend(cursor_t const& cur) const noexcept {
        if (dir == dir_t::Right) return cur.extendRight();
        return cur.extendLeft();
    }

    auto extendKStep(cursor_t const& cur) const noexcept {
        if constexpr (HasKStep) {
            if (dir == dir_t::Right) return cur.extendRightKStep();
            return cur.extendLeftKStep();
        } else {
            return extend(cur);
        }
    }


    auto computeErrorCases() const {
        auto const tinfo               = side->info;
        bool const mismatchAllowed     = e+1 <= ub;
        if (!mismatchAllowed) return std::make_tuple(false, false, false, false);
        bool const substitutionAllowed = lb <= e+partitionPart and e+1 <= ub;
        bool const insertionAllowed    = Edit && tinfo != 'S' && tinfo != 'D' && substitutionAllowed;
        bool const deletionAllowed     = Edit && tinfo != 'S' && tinfo != 'I' && mismatchAllowed;
        return std::make_tuple(mismatchAllowed, substitutionAllowed, insertionAllowed, deletionAllowed);
    }

    auto computeMatchCase(size_t nextSymb) const {
        auto const tinfo        = side->info;
        bool const matchAllowed = (partitionPart > 1 or lb <= e)
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
            return clearCursor(cur, [&](auto const& cur) {
                auto il = sideImpl[0].info;
                auto ir = sideImpl[2].info;
                if (!Edit || ((il == 'M' or il == 'I') and (ir == 'M' or ir == 'I'))) {
                    if (search.l.back() <= e and e <= search.u.back()) {
                        return delegate(cur, e);
                    }
                }
                return false;
            });
        }

        bool right = (part==0) || (search.pi[part-1] < search.pi[part]);
        auto r_side = Restore{side, right?&sideImpl[2]:&sideImpl[0]};
        auto r_part = Restore{partitionPart, partition[search.pi[part]]};
        auto r_lb   = Restore{lb, search.l[part]};
        auto r_ub   = Restore{ub, search.u[part]};

        auto newDir = right?dir_t::Right:dir_t::Left;
        if (newDir == dir) {
            return search_next_symb(cur);
        } else {
            // clearCursor and search other direction
            return clearCursor(cur, [&](auto const& _cur) {
                auto r_dir = Restore{dir, newDir};
                return search_next_symb(_cur);
            });
        }
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

    template <bool ExtendCursor = true>
    bool search_next_symb(cursor_t const& cur) const {
        if (cur.count() == 1) {
            return clearCursor(cur, [&](auto const& cur) {
                return search_next_symb_single(cur);
            });
        }

        if (side->usedSymbBuffer % KStep == 0) {
            auto nextSymb = query[side->queryPos];

            auto matchAllowed = computeMatchCase(nextSymb);
            auto [mismatchAllowed, substitutionAllowed, insertionAllowed, deletionAllowed] = computeErrorCases();

            if (!mismatchAllowed && matchAllowed) {
                return search_next_part_no_errors(cur);
            }
        }
        // Is the next step reached?
        if constexpr (ExtendCursor) {
            if (side->usedSymbBuffer % KStep == 0) {
                // require expansion of cursor buffer
                auto cursors = extendKStep(cur);
                auto r_nextCursor = Restore{nextCursors, &cursors};
                return search_next_symb<false>(cur);
            }
        }


        // will not need a new extension
        if (side->usedSymbBuffer % KStep + 1 < KStep) {
            auto r_i    = Restore{side->info};
            auto r_lqr  = Restore{side->lastQRank};
            auto r_qp   = Restore{side->queryPos};
            auto r_p    = Restore{partitionPart};
            auto r_e    = Restore{e};
            auto r_used = Restore{side->usedSymbBuffer};

            side->usedSymbBuffer += 1;
            while (true) {
                auto nextSymb = query[side->queryPos];

                auto matchAllowed = computeMatchCase(nextSymb);
                auto [mismatchAllowed, substitutionAllowed, insertionAllowed, deletionAllowed] = computeErrorCases();

                if (!mismatchAllowed && matchAllowed) {
                    auto r_lr  = Restore{side->lastRank, nextSymb};
                    auto r_lqr = Restore{side->lastQRank, nextSymb};
                    side->info = 'M';
                    side->symbBuffer[side->usedSymbBuffer-1] = nextSymb;

                    // copied from check_and_search_next_symb
                    auto r_qp = RestoreAdd{side->queryPos, dir};
                    auto r_p  = RestoreSub{partitionPart, 1};

                    if (partitionPart == 0) {
                        return search_next_part(cur);
                    }

                    return search_next_part_no_errors(cur);
                }


                if (matchAllowed) {
                    auto r_lr  = Restore{side->lastRank, nextSymb};
                    auto r_lqr = Restore{side->lastQRank, nextSymb};
                    side->info = 'M';
                    side->symbBuffer[side->usedSymbBuffer-1] = nextSymb;
                    auto f = check_and_search_next_symb(cur);
                    if (f) return true;
                }
                if (!mismatchAllowed) break;
                e += 1;
                for (uint64_t i{FirstSymb}; i < Sigma; ++i) {
                    side->symbBuffer[side->usedSymbBuffer-1] = i;
                    auto r_lr = Restore{side->lastRank, i};
                    if (deletionAllowed) {
                        side->info = 'D';
                        auto f = search_next_symb(cur); // deletion occurred in query
                        if (f) return true;
                    }
                    if (!substitutionAllowed) continue;
                    if (i == nextSymb) continue;

                    auto r_lqr = Restore{side->lastQRank, nextSymb};
                    side->info = 'S';
                    auto f = check_and_search_next_symb(cur);
                    if (f) return true;
                }

                if (!insertionAllowed) break;
                side->lastQRank = nextSymb;
                side->info = 'I';
                side->queryPos += dir;
                partitionPart -= 1;
                if (partitionPart == 0) {
                    side->usedSymbBuffer -= 1;
                    return search_next_part(cur);
                }
            }
            return false;
        } else { // new extension will be required
            auto r_i    = Restore{side->info};
            auto r_lqr  = Restore{side->lastQRank};
            auto r_qp   = Restore{side->queryPos};
            auto r_p    = Restore{partitionPart};
            auto r_e    = Restore{e};
            auto r_used = Restore{side->usedSymbBuffer};

            auto const& cursors = *nextCursors;

            side->usedSymbBuffer += 1;
            auto kstepSymb = size_t{};
            for (size_t i{0}; i < KStep-1; ++i) {
                kstepSymb = kstepSymb*Sigma + side->symbBuffer[side->usedSymbBuffer-KStep+i];
            }
            kstepSymb = kstepSymb*Sigma;
            while (true) {
                auto nextSymb = query[side->queryPos];

                auto matchAllowed = computeMatchCase(nextSymb);
                auto [mismatchAllowed, substitutionAllowed, insertionAllowed, deletionAllowed] = computeErrorCases();

                if (!mismatchAllowed && matchAllowed) {
                    auto const& newCur = cursors[kstepSymb + nextSymb];
                    if (newCur.count() == 0) return false;

                    auto r_lr  = Restore{side->lastRank, nextSymb};
                    auto r_lqr = Restore{side->lastQRank, nextSymb};
                    side->info = 'M';

                    // copied from check_and_search_next_symb
                    auto r_qp = RestoreAdd{side->queryPos, dir};
                    auto r_p  = RestoreSub{partitionPart, 1};

                    if (partitionPart == 0) {
                        return search_next_part(newCur);
                    }

                    return search_next_part_no_errors(newCur);
                }

                if (matchAllowed) {
                    auto const& newCur = cursors[kstepSymb + nextSymb];

                    auto r_lr  = Restore{side->lastRank, nextSymb};
                    auto r_lqr = Restore{side->lastQRank, nextSymb};
                    side->info = 'M';
                    auto f = check_and_search_next_symb(newCur);
                    if (f) return true;
                }
                if (!mismatchAllowed) break;
                e += 1;
                for (uint64_t i{FirstSymb}; i < Sigma; ++i) {
                    auto const& newCur = cursors[kstepSymb + i];
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
                    side->usedSymbBuffer -= 1;
                    return search_next_part(cur);
                }
            }
            return false;
        }
    }

    /** Searches the (rest of) a part without any errors
     *
     * May only be called if no errors are allowed
     */
    bool search_next_part_no_errors(cursor_t const& _cur) const {
        auto r_qp   = Restore{side->queryPos};
        assert(e == ub);

        return clearCursor(_cur, [&](auto cur) {
            assert(side->usedSymbBuffer % KStep == 0);

            auto kstep_loops = partitionPart / KStep;               // number of k-steps
            auto loop_tail   = partitionPart - kstep_loops * KStep; // trailing number of single steps

            auto buffer = std::array<size_t, KStep>{};
            for (size_t i{0}; i < kstep_loops; ++i) {
                for (size_t j{0}; j < KStep; ++j) {
                    buffer[j] = query[side->queryPos];
                    side->queryPos = side->queryPos + dir;
                }
                cur = extendKStep(cur, buffer);
                if (cur.count() == 0) return false; // early abort, if results are empty
            }

            auto nextSymb = buffer.back();
            for (size_t i{0}; i < loop_tail; ++i) {
                nextSymb = query[side->queryPos];
                side->queryPos = side->queryPos + dir;
                cur = extend(cur, nextSymb);
                if (cur.count() == 0) return false; // early abort, if results are empty
            }
            auto r_lr   = Restore{side->lastRank, nextSymb};
            auto r_lqr  = Restore{side->lastQRank, nextSymb};
            auto r_p    = Restore{partitionPart, 0};
            auto r_i    = Restore{side->info, 'M'};
            return search_next_part(cur);
        });
    }

    bool search_next_symb_single(cursor_t const& cur) const {
        auto r_e    = Restore{e};
        auto r_lqr  = Restore{side->lastQRank};
        auto r_i    = Restore{side->info};
        auto r_qp   = Restore{side->queryPos};
        auto r_p    = Restore{partitionPart};

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
