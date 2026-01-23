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

template <bool Edit, typename index_t, typename query_t, typename search_t, typename report_t>
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
    report_t const& report;

    enum dir_t : uint8_t {
        Left =  0,
        Right = 1,
    };

    struct Node {
        struct Side {
            uint8_t lastIRank{};
            uint8_t lastQRank{};
            char    info{'M'};
            size_t  queryPos{};

        };
        std::array<Side, 2> side;
        size_t e{};
        size_t part{std::numeric_limits<size_t>::max()}; // first action will be overflow, this is desired effect
        dir_t dir{dir_t::Right};
        size_t partitionPart{};
    };

    Node rootNode{};

    Search(index_t const& _index, query_t const& _query, search_t const& _search, std::vector<size_t> const& _partition, report_t const& _report)
        : index     {_index}
        , query     {_query}
        , search    {_search}
        , partition {_partition}
        , report  {_report}
    {
        // check how many characters are before the first query char
        for (size_t i{0}; i < search.pi[0]; ++i) {
            rootNode.side[0].queryPos += partition[i];
            rootNode.side[1].queryPos += partition[i];
        }

        // Notice: rootNode.side[0].queryPos might 'underflow' (also at other spots)
        // but that is fine, after an underflow this variable is not
        // being accessed before overflowing again.
        rootNode.side[0].queryPos -= 1;
    }

    auto computeErrorCases(Node const& node) const {
        auto lb = search.l[node.part];
        auto ub = search.u[node.part];
        auto const tinfo               = node.side[node.dir].info;
        bool const mismatchAllowed     = node.e+1 <= ub;
        if (!mismatchAllowed) return std::make_tuple(false, false, false, false);
        bool const substitutionAllowed = ((Edit && (node.partitionPart > 0 || lb <= node.e+1)) || (!Edit && lb < node.e+node.partitionPart)) && mismatchAllowed;
        bool const insertionAllowed    = Edit && tinfo != 'S' && tinfo != 'D' && substitutionAllowed;
        bool const deletionAllowed     = Edit && tinfo != 'S' && tinfo != 'I' && mismatchAllowed;
        return std::make_tuple(mismatchAllowed, substitutionAllowed, insertionAllowed, deletionAllowed);
    }

    auto computeMatchCase(Node const& node, size_t nextSymb) const {
        auto lb = search.l[node.part];
        auto ub = search.u[node.part];

        auto const tinfo        = node.side[node.dir].info;
        bool const matchAllowed = (!Edit && (node.partitionPart > 0 || lb <= node.e)
                                    || (Edit && lb <= node.e + node.partitionPart))
                                 and node.e <= ub
                                 and (tinfo != 'I' or nextSymb != node.side[node.dir].lastQRank)
                                 and (tinfo != 'D' or nextSymb != node.side[node.dir].lastIRank);
        return matchAllowed;
    }

    auto computeNoErrors(Node const& node, size_t nextSymb) const {
        auto matchAllowed = computeMatchCase(node, nextSymb);
        auto ub = search.u[node.part];
        return node.e == ub && matchAllowed;
    }

    auto extend(cursor_t const& cur, uint64_t symb, dir_t dir) const noexcept {
        if (dir == dir_t::Right) {
            return cur.extendRight(symb);
        } else {
            return cur.extendLeft(symb);
        }
    }

/*    auto extendKStep(cursor_t const& cur, std::span<size_t const, KStep> symb) const noexcept {
        if constexpr (HasKStep) {
            if (dir == dir_t::Right) return cur.extendRightKStep(symb);
            else return cur.extendLeftKStep(symb);
        } else {
            static_assert(KStep == 1);
            return extend(cur, symb[0]);
        }
    }*/

    auto extend(cursor_t const& cur, dir_t dir) const noexcept {
        if (dir == dir_t::Right) return cur.extendRight();
        return cur.extendLeft();
    }

/*    auto extendKStep(cursor_t const& cur) const noexcept {
        if constexpr (HasKStep) {
            if (dir == dir_t::Right) return cur.extendRightKStep();
            return cur.extendLeftKStep();
        } else {
            return extend(cur);
        }
    }*/

    std::vector<Node> stackNode;
    bool run() {
        stackNode.clear();
        auto This = this; // required by gcc 15.2
        auto handleNode = [&](this auto const& self, cursor_t const& cur, Node node) -> bool {
            // check if end of part is reached, move to next part if required
            if (node.partitionPart == 0) {
                node.part += 1;

                auto const& search    = This->search;
                auto const& partition = This->partition;

                // abort if end of all partitions is reached
                if (node.part == partition.size()) {
                    auto il = node.side[0].info;
                    auto ir = node.side[1].info;
                    if (!Edit || ((il == 'M' or il == 'I') and (ir == 'M' or ir == 'I'))) {
                        if (search.l.back() <= node.e and node.e <= search.u.back()) {
                            return This->report(cur, node.e);
                        }
                    }
                    return false;
                }
                auto rightDir = (node.part == 0) || (search.pi[node.part-1] < search.pi[node.part]);
                node.partitionPart = partition[search.pi[node.part]];
                node.dir = rightDir?dir_t::Right:dir_t::Left;
            }
            auto dirStep = node.dir == dir_t::Right?1:-1;

            // Optimization: if no errors are allowed, extend cursor directly
            auto nextSymb = This->query[node.side[node.dir].queryPos];
            auto noErrors = This->computeNoErrors(node, nextSymb);
            if (noErrors) {
                auto& side = node.side[node.dir];
                auto newCur = cur;
                auto loops = node.partitionPart;
                auto dir = node.dir;
                for (size_t i{0}; i < loops; ++i) {
                    auto pos = side.queryPos + i*dirStep;
                    nextSymb = This->query[pos];
                    if (newCur.count() == 1) {
                        if (dir == dir_t::Right) { if (newCur.symbolRight() != nextSymb) return false; }
                        else                     { if (newCur.symbolLeft() != nextSymb) return false; }
                        if constexpr (requires() { { newCur.extendRightBySymbol() }; }) {
                            if (dir == dir_t::Right) newCur = newCur.extendRightBySymbol(nextSymb);
                            else                     newCur = newCur.extendLeftBySymbol(nextSymb);
                        } else {
                            if (dir == dir_t::Right) newCur = newCur.extendRight(nextSymb);
                            else                     newCur = newCur.extendLeft(nextSymb);
                        }
                    } else {
                        newCur = This->extend(newCur, nextSymb, dir);
                    }
                    if (newCur.count() == 0) return false; // early abort, if results are empty
                }
                node.partitionPart = 0;
                side.queryPos += loops * dirStep;
                side.info = 'M';
                side.lastQRank = nextSymb;
                side.lastIRank = nextSymb;
                return self(newCur, node);
            }

            // Optimization: if cursor only covers a single row, we only need to expand that
            if (cur.count() == 1) {
                // Optimization: if cursor only covers a single row, we only need to expand that
                auto [curISymb, icursorNext] = [&]() -> std::tuple<size_t, cursor_t> {
                    if constexpr (requires() { { cur.extendRightBySymbol() }; }) {
                        if (node.dir == dir_t::Right) {
                            return cur.extendRightBySymbol();
                        } else {
                            return cur.extendLeftBySymbol();
                        }
                    }
                    if (node.dir == dir_t::Right) {
                        auto symb = cur.symbolRight();
                        auto cur_ = cur.extendRight(symb);
                        return {symb, cur_};
                    } else {
                        auto symb = cur.symbolLeft();
                        auto cur_ = cur.extendLeft(symb);
                        return {symb, cur_};
                    }
                }();
                auto optInsert = std::optional<Node>{};
                bool abort{};
                This->expandNode(node, [&](Node const& child) {
                    if (abort) return;
                    if (child.side[node.dir].info == 'I') {
                        optInsert = child;
                        return;
                    }

                    if (child.side[node.dir].info == 'M') {
                        if (child.side[node.dir].lastIRank != curISymb) return;
                        abort = self(icursorNext, child);
                    } else if (child.side[node.dir].info == 'S') {
                        if (curISymb == child.side[node.dir].lastQRank) return;
                        auto _child = child;
                        _child.side[node.dir].lastIRank = curISymb;
                        abort = self(icursorNext, _child);
                    } else if (child.side[node.dir].info == 'D') {
                        auto _child = child;
                        _child.side[node.dir].lastIRank = curISymb;
                        abort = self(icursorNext, _child);
                    }
                });
                if (abort) return true;
                if (optInsert) {
                    auto const& child = *optInsert;
                    // inserted into query
                    if (auto f = self(cur, child)) return true;
                }
                return false;
            }


            // standard, expand node by following the direct connected edges
            auto startIdx = This->stackNode.size();
            auto optInsert = std::optional<Node>{};
            bool errorMatch{};
            This->expandNode(node, [&](Node const& child) {
                if (child.side[node.dir].info == 'I') {
                    optInsert = child;
                } else {
                    if (child.side[node.dir].info != 'M') errorMatch = true;
                    This->stackNode.push_back(child);
                }
            });
            auto endIdx = This->stackNode.size();
            assert(endIdx >= startIdx);

            if (!errorMatch) { //!TODO magic constant of 2, probably should be dependent of Sigma
                // Optimization: if only 1-2 paths are followed
                // do not expand cursor to all possibilities
                for (size_t i{startIdx}; i < endIdx; ++i) {
                    auto const& child = This->stackNode[i];

                    // match/substitution
                    auto newCur = This->extend(cur, child.side[node.dir].lastIRank, node.dir);
                    if (newCur.count() == 0) continue;
                    if (auto f = self(newCur, child)) return true;
                }
            } else {
                // Expand cursor to all paths
                auto cursors = This->extend(cur, node.dir);
                for (size_t i{startIdx}; i < endIdx; ++i) {
                    auto& child = This->stackNode[i];

                    if (child.side[node.dir].info == 'M') {
                        auto const& newCur = cursors[child.side[node.dir].lastIRank];
                        if (newCur.count() == 0) continue;
                        if (auto f = self(newCur, child)) return true;
                    } else if (child.side[node.dir].info == 'S') {
                        for (uint64_t i{FirstSymb}; i < Sigma; ++i) {
                            if (i == child.side[node.dir].lastQRank) continue;
                            auto const& newCur = cursors[i];
                            if (newCur.count() == 0) continue;
                            child.side[node.dir].lastIRank = i;
                            if (auto f = self(newCur, child)) return true;
                        }
                    } else if (child.side[node.dir].info == 'D') {
                        for (uint64_t i{FirstSymb}; i < Sigma; ++i) {
                            auto const& newCur = cursors[i];
                            if (newCur.count() == 0) continue;
                            child.side[node.dir].lastIRank = i;
                            if (auto f = self(newCur, child)) return true;
                        }
                    }
                }
            }

            if (optInsert) {
                auto const& child = *optInsert;
                // inserted into query
                if (auto f = self(cur, child)) return true;
            }
            This->stackNode.resize(startIdx);
            return false;
        };
        auto cur = cursor_t{index};
        return handleNode(cur, rootNode);
    }
    void expandNode(Node node, auto const& cb) {
        auto dirStep = node.dir == dir_t::Right?1:-1;
        auto* side = &node.side[node.dir];
        auto nextSymb = query[side->queryPos];
        side->queryPos += dirStep;
        node.partitionPart -= 1;

        while (true) {
            auto matchAllowed = computeMatchCase(node, nextSymb);
            auto [mismatchAllowed, substitutionAllowed, insertionAllowed, deletionAllowed] = computeErrorCases(node);

            if (matchAllowed) {
                auto r_lr  = Restore{side->lastIRank, nextSymb};
                auto r_lqr = Restore{side->lastQRank, nextSymb};
                side->info = 'M';
                cb(node);
            }
            if (!mismatchAllowed) return;
            node.e += 1;

            if (deletionAllowed) {
                auto r_pos   = RestoreSub{side->queryPos, dirStep};
                auto r_part  = RestoreAdd{node.partitionPart, 1};
                side->info = 'D';
                cb(node); // deletion occurred in query
            }
            if (substitutionAllowed) {
                auto r_lqr   = Restore{side->lastQRank, nextSymb};
                side->info = 'S';
                cb(node); // character in query was substituted
            }

            if (!insertionAllowed) return;
            side->lastQRank = nextSymb;
            side->info = 'I';

            if (node.partitionPart == 0) {
                // insertion does not trigger a direct expansion of the node
                cb(node);
                return;
            }

            nextSymb = query[side->queryPos];
            side->queryPos += dirStep;
            node.partitionPart -= 1;
        }
    }
};


template <bool Edit, typename index_t, Sequence query_t, typename report_t>
void search_impl(index_t const& index, query_t const& query, search_scheme::Scheme const& search_scheme, std::vector<size_t> const& partition, report_t&& report) {
    using cursor_t = select_cursor_t<index_t>;
    using R = std::decay_t<decltype(report(std::declval<cursor_t>(), 0))>;

    auto internal_report = [&]() {
        if constexpr (std::same_as<R, bool>) {
            return report;
        } else {
            return [=](auto cur, size_t e) {
                report(cur, e);
                return std::false_type{};
            };
        }
    }();

    for (auto const& search : search_scheme) {
        bool f = Search<Edit, index_t, query_t, decltype(search), decltype(internal_report)>{index, query, search, partition, internal_report}.run();
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
 * \param report_t: callback function to report the results, Must accept there parameters: size_t qidx, auto cur, size_t e
 *          qidx: the position of the query inside the queries range
 *          cur:  cursor of the FMIndex marking the result positions
 *          e:    number of errors that these matches have
 * \param searchSelectSearchScheme_t: callback that helps selecting a proper search scheme, Must accept one parameters: size_t length
 *          length: length of query
 */
template <bool Edit, typename index_t, Sequences queries_t, typename selectSearchScheme_t, typename report_t>
void search_n_impl(index_t const& index, queries_t&& queries, selectSearchScheme_t&& selectSearchScheme, std::vector<size_t> const& partition, report_t&& report, size_t n) {
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
            report(qidx, cur, e);
            return ct == n;
        });
    }
}


// convenience function, with passed search scheme and multiple queries
template <bool Edit=true, typename index_t, Sequences queries_t, typename report_t>
void search(index_t const& index, queries_t&& queries, search_scheme::Scheme const& search_scheme, std::vector<size_t> const& partition, report_t&& report, size_t n = std::numeric_limits<size_t>::max()) {
    // function that selects a search scheme
    auto selectSearchScheme = [&]([[maybe_unused]] size_t length) -> auto& {
        return search_scheme;
    };
    search_n_impl<Edit>(index, queries, selectSearchScheme, partition, report, n);
}

}
