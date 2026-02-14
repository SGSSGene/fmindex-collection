// SPDX-FileCopyrightText: 2026 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "CachedSearchScheme.h"
#include "Restore.h"
#include "SelectCursor.h"

#include <array>
#include <cstddef>

//#define opt1 1
//#define opt1_1r 0
//#define opt2 0
//#define opt3 0 // not available for == 1
/**
 * like search_ng28 but:
 *  - optimized for KStep approach
 */
namespace fmc::search_ng28_kstep {

struct Optimizations {
    size_t opt1{2};
    bool opt1_1r{true};
    size_t opt2{2};
    size_t opt3{2};
};

template <bool Edit, typename index_t, typename query_t, typename search_t, typename report_t, Optimizations OPT>
struct Search {
    static_assert(OPT.opt1 == 0 || OPT.opt1 == 1 || OPT.opt1 == 2);
    static_assert(OPT.opt2 == 0 || OPT.opt2 == 1 || OPT.opt2 == 2);
    static_assert(OPT.opt3 == 0 || OPT.opt3 == 2);

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
        bool const substitutionAllowed = (((Edit && (node.partitionPart > 0 || lb <= node.e+1)) || (!Edit && lb < node.e+node.partitionPart)));
        bool const insertionAllowed    = Edit && tinfo != 'S' && tinfo != 'D' && substitutionAllowed;
        bool const deletionAllowed     = Edit && tinfo != 'S' && tinfo != 'I';
        return std::make_tuple(mismatchAllowed, substitutionAllowed, insertionAllowed, deletionAllowed);
    }

    auto computeMatchCase(Node const& node, size_t nextSymb) const {
        auto lb = search.l[node.part];
        auto ub = search.u[node.part];

        auto const tinfo        = node.side[node.dir].info;
        bool const matchAllowed = ((!Edit && (node.partitionPart > 0 || lb <= node.e))
                                    || (Edit && lb <= node.e + node.partitionPart))
                                 and node.e <= ub
                                 and (tinfo != 'I' or nextSymb != node.side[node.dir].lastQRank)
                                 and (tinfo != 'D' or nextSymb != node.side[node.dir].lastIRank);
        return matchAllowed;
    }

    auto computeMatchCaseSimple(Node const& node) const {
        auto lb = search.l[node.part];
        auto ub = search.u[node.part];

        auto const tinfo        = node.side[node.dir].info;
        bool const matchAllowed = ((!Edit && (node.partitionPart > 0 || lb <= node.e))
                                    || (Edit && lb <= node.e + node.partitionPart))
                                 and node.e <= ub;
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

    auto extendKStep(cursor_t const& cur, std::span<size_t const, KStep> symb, dir_t dir) const noexcept requires HasKStep {
        if (dir == dir_t::Right) return cur.extendRightKStep(symb);
        else return cur.extendLeftKStep(symb);
    }

    auto extend(cursor_t const& cur, dir_t dir) const noexcept {
        if (dir == dir_t::Right) return cur.extendRight();
        return cur.extendLeft();
    }

    auto extendKStep(cursor_t const& cur, dir_t dir) const noexcept requires HasKStep {
        if (dir == dir_t::Right) return cur.extendRightKStep();
        return cur.extendLeftKStep();
    }

    bool run() {
        if constexpr (HasKStep) {
            return run_2step();
        } else {
            return run_1step();
        }
    }

    std::vector<Node> stackNode;
    bool run_1step() {
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
            if constexpr (OPT.opt1 > 0) {
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

                        // Optimization: if no errors and only a single row is available
                        if constexpr (OPT.opt1_1r) {
                            if (newCur.count() == 1) {
                                auto s = (dir == dir_t::Right)?newCur.symbolRight():newCur.symbolLeft();
                                if (s != nextSymb) return false;
                                if constexpr (requires() { { newCur.extendRightBySymbol(nextSymb) }; }) {
                                    newCur = (dir == dir_t::Right)?newCur.extendRightBySymbol(nextSymb):newCur.extendLeftBySymbol(nextSymb);
                                } else {
                                    newCur = (dir == dir_t::Right)?newCur.extendRight(nextSymb):newCur.extendLeft(nextSymb);
                                }
                            } else {
                                newCur = This->extend(newCur, nextSymb, dir);
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
            }

            if constexpr (OPT.opt2 > 0) {
                // Optimization: if cursor only covers a single row, we only need to expand that
                if (cur.count() == 1) {
                    // Optimization: if cursor only covers a single row, we only need to expand that
                    auto [curISymb, icursorNext] = [&]() -> std::tuple<size_t, cursor_t> {
                        if constexpr (requires() { { cur.extendRightBySymbol() }; }) {
                            return (node.dir == dir_t::Right)?cur.extendRightBySymbol():cur.extendLeftBySymbol();
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
                    This->expandNode(node, [&](Node& child) {
                        if (abort) return;
                        switch(child.side[node.dir].info) {
                        case 'I':
                            optInsert = child;
                            return;
                        case 'M':
                            if (child.side[node.dir].lastIRank != curISymb) return;
                            abort = self(icursorNext, child);
                            break;
                        case 'S':
                            if (curISymb == child.side[node.dir].lastQRank) return;
                        case 'D': // fallthrough
                        default: // not reachable
                            {
                                auto r_irank = Restore{child.side[node.dir].lastIRank, curISymb};
                                abort = self(icursorNext, child);
                            }
                            break;
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

            if (!errorMatch) {
                // Optimization: if only matches exist
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
                    switch(child.side[node.dir].info) {
                    case 'M': {
                        auto const& newCur = cursors[child.side[node.dir].lastIRank];
                        if (newCur.count() == 0) continue;
                        if (auto f = self(newCur, child)) return true;
                        }
                        break;
                    case 'S':
                        for (uint64_t i{FirstSymb}; i < Sigma; ++i) {
                            if (i == child.side[node.dir].lastQRank) continue;
                            auto const& newCur = cursors[i];
                            if (newCur.count() == 0) continue;
                            child.side[node.dir].lastIRank = i;
                            if (auto f = self(newCur, child)) return true;
                        }
                        break;
                    case 'D':
                        for (uint64_t i{FirstSymb}; i < Sigma; ++i) {
                            auto const& newCur = cursors[i];
                            if (newCur.count() == 0) continue;
                            child.side[node.dir].lastIRank = i;
                            if (auto f = self(newCur, child)) return true;
                        }
                        break;
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

    /**
     * If SimpleError is set to true, the previous character for deciding if a match is allowed, is being ignored
     */
    template <bool SimpleError = false>
    void expandNode(Node node, auto const& cb) {
        auto dirStep = node.dir == dir_t::Right?1:-1;
        auto* side = &node.side[node.dir];
        auto nextSymb = query[side->queryPos];
        side->queryPos += dirStep;
        node.partitionPart -= 1;

        while (true) {
            auto matchAllowed = SimpleError?computeMatchCaseSimple(node):computeMatchCase(node, nextSymb);
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

    std::vector<std::tuple<Node, Node>> stackNode2; // buffer for symbols of each node
    bool run_2step() requires HasKStep {
        stackNode.clear();
        stackNode2.clear();
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

            if constexpr (OPT.opt1 == 1) {
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

                        // Optimization: if no errors and only a single row is available
                        if constexpr (OPT.opt1_1r) {
                            if (newCur.count() == 1) {
                                auto s = (dir == dir_t::Right)?newCur.symbolRight():newCur.symbolLeft();
                                if (s != nextSymb) return false;
                                if constexpr (requires() { { newCur.extendRightBySymbol(nextSymb) }; }) {
                                    newCur = (dir == dir_t::Right)?newCur.extendRightBySymbol(nextSymb):newCur.extendLeftBySymbol(nextSymb);
                                } else {
                                    newCur = (dir == dir_t::Right)?newCur.extendRight(nextSymb):newCur.extendLeft(nextSymb);
                                }
                            } else {
                                newCur = This->extend(newCur, nextSymb, dir);
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
            }
            if constexpr (OPT.opt1 == 2) {
                // Optimization: if no errors are allowed, extend cursor directly
                auto nextSymb = This->query[node.side[node.dir].queryPos];
                auto noErrors = This->computeNoErrors(node, nextSymb);
                if (noErrors) {
                    auto& side = node.side[node.dir];
                    auto newCur = cur;
                    auto loops = node.partitionPart/KStep;
                    auto dir = node.dir;
                    auto offset = node.partitionPart % KStep;

                    for (size_t i{0}; i < loops; ++i) {
                        auto pos1 = side.queryPos + (i*KStep)*dirStep;
                        auto pos2 = side.queryPos + (i*KStep+1)*dirStep;
                        auto s1 = This->query[pos1];
                        auto s2 = This->query[pos2];

                        //!TODO only works for KStep==2
                        auto symbs = std::array<size_t, KStep> {
                            s1,
                            s2
                        };
                        auto kSymb = s1 * Sigma + s2;


                        if constexpr (OPT.opt1_1r) {
                            // Optimization if only a single row exists
                            if (newCur.count() == 1) {
                                auto s = (dir == dir_t::Right)?newCur.symbolRightKStep():newCur.symbolLeftKStep();
                                if (s != kSymb) return false;
                                if constexpr (requires() { { newCur.extendRightBySymbolKStep(kSymb) }; }) {
                                    newCur = (dir == dir_t::Right)?newCur.extendRightBySymbolKStep(kSymb):newCur.extendLeftBySymbolKStep(kSymb);
                                } else {
                                    newCur = (dir == dir_t::Right)?newCur.extendRightKStep(symbs):newCur.extendLeftKStep(symbs);
                                }
                            } else {
                                newCur = This->extendKStep(newCur, symbs, dir);
                            }
                        } else {
                            newCur = This->extendKStep(newCur, symbs, dir);
                        }

                        if (newCur.count() == 0) return false; // early abort, if results are empty
                        nextSymb = s2;
                    }
                    side.queryPos += (loops * KStep) * dirStep;
                    for (size_t i{0}; i < offset; ++i) {
                        if (newCur.count() == 0) return false;
                        nextSymb = This->query[side.queryPos];
                        newCur = This->extend(newCur, nextSymb, dir);
                        side.queryPos += dirStep;
                    }

                    node.partitionPart = 0;
                    side.info = 'M';
                    side.lastQRank = nextSymb;
                    side.lastIRank = nextSymb;
                    return self(newCur, node);
                }
            }

            if constexpr (OPT.opt2 == 1) {
                // Optimization: if cursor only covers a single row, we only need to expand that
                if (cur.count() == 1) {
                    // Optimization: if cursor only covers a single row, we only need to expand that
                    auto [curISymb, icursorNext] = [&]() -> std::tuple<size_t, cursor_t> {
                        if constexpr (requires() { { cur.extendRightBySymbol() }; }) {
                            return (node.dir == dir_t::Right)?cur.extendRightBySymbol():cur.extendLeftBySymbol();
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
                    This->expandNode(node, [&](Node& child) {
                        if (abort) return;
                        switch(child.side[node.dir].info) {
                        case 'I':
                            optInsert = child;
                            return;
                        case 'M':
                            if (child.side[node.dir].lastIRank != curISymb) return;
                            abort = self(icursorNext, child);
                            break;
                        case 'S':
                            if (curISymb == child.side[node.dir].lastQRank) return;
                        case 'D': // fallthrough
                        default: // not reachable
                            {
                                auto r_irank = Restore{child.side[node.dir].lastIRank, curISymb};
                                abort = self(icursorNext, child);
                            }
                            break;
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
            }

            if constexpr (OPT.opt2 == 2) {
                // Optimization: if cursor only covers a single row, we only need to expand that
                if (cur.count() == 1) {
                    // Optimization: if cursor only covers a single row, we only need to expand that
                    auto [curISymb, icursorNext] = [&]() -> std::tuple<size_t, cursor_t> {
                        if constexpr (requires() { { cur.extendRightBySymbolKStep() }; }) {
                            return (node.dir == dir_t::Right)?cur.extendRightBySymbolKStep():cur.extendLeftBySymbolKStep();
                        }
                        if (node.dir == dir_t::Right) {
                            auto symb = cur.symbolRightKStep();
                            auto cur_ = cur.extendRightKStep(symb);
                            return {symb, cur_};
                        } else {
                            auto symb = cur.symbolLeftKStep();
                            auto cur_ = cur.extendLeftKStep(symb);
                            return {symb, cur_};
                        }
                    }();
                    auto optInsert = std::optional<Node>{};
                    bool abort{};
                    This->expandNode(node, [&](Node& child) {
                        if (abort) return;
                        switch(child.side[node.dir].info) {
                        case 'I':
                            optInsert = child;
                            return;
                        case 'M':
                            if (child.side[node.dir].lastIRank != curISymb) return;
                            abort = self(icursorNext, child);
                            break;
                        case 'S':
                            if (curISymb == child.side[node.dir].lastQRank) return;
                        case 'D': // fallthrough
                        default: // not reachable
                            {
                                auto r_irank = Restore{child.side[node.dir].lastIRank, curISymb};
                                abort = self(icursorNext, child);
                            }
                            break;
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
            }

            if constexpr (OPT.opt3 == 2) {
                // standard, expand node by following the direct connected edges
                auto startIdx = This->stackNode.size();
                auto startIdx2 = This->stackNode2.size();

    //            auto r_stack  = RestoreCB{[=]() { This->stackNode.resize(startIdx); }};
    //            auto r_stack2 = RestoreCB{[=]() { This->stackNode2.resize(startIdx2); }};
                //This->stackNode.resize(startIdx);
                //This->stackNode2.resize(startIdx2);




                std::vector<Node> optInsert;
                std::vector<std::tuple<Node, Node>> optInsert2;
    //            auto optInsert = std::optional<Node>{};
                bool errorMatch{};
                bool errorMatch2{};
                This->expandNode_2Step(node, [&](Node const& child) {
                    if (child.side[node.dir].info == 'I') {
                        optInsert.push_back(child);
    //                    optInsert = child;
                        return;
                    }
                    if (child.side[node.dir].info != 'M') errorMatch = true;
                    This->stackNode.push_back(child);
                }, [&](Node const& child1, Node const& child2) {
                    if (child1.side[node.dir].info == 'I') {
                        throw std::runtime_error{"this should not happen 4"};
                    }
                    if (child2.side[node.dir].info == 'I') {
                        optInsert2.emplace_back(child1, child2);
                        //!TODO optInsert2 is not handled anywhere
                        return;
                    }
                    if (child1.side[node.dir].info != 'M'
                        || child2.side[node.dir].info != 'M') {
                        errorMatch2 = true;
                    }
                    This->stackNode2.emplace_back(child1, child2);
                });
                auto endIdx = This->stackNode.size();
                auto endIdx2 = This->stackNode2.size();
                assert(endIdx >= startIdx);
                assert(endIdx2 >= startIdx2);

               if (endIdx > startIdx) {
                    if (!errorMatch) {
                        // Optimization: if only matches exist
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
                            switch(child.side[node.dir].info) {
                            case 'M': {
                                auto const& newCur = cursors[child.side[node.dir].lastIRank];
                                if (newCur.count() == 0) continue;
                                if (auto f = self(newCur, child)) return true;
                                }
                                break;
                            case 'S':
                                for (uint64_t i{FirstSymb}; i < Sigma; ++i) {
                                    if (i == child.side[node.dir].lastQRank) continue;
                                    auto const& newCur = cursors[i];
                                    if (newCur.count() == 0) continue;
                                    child.side[node.dir].lastIRank = i;
                                    if (auto f = self(newCur, child)) return true;
                                }
                                break;
                            case 'D':
                                for (uint64_t i{FirstSymb}; i < Sigma; ++i) {
                                    auto const& newCur = cursors[i];
                                    if (newCur.count() == 0) continue;
                                    child.side[node.dir].lastIRank = i;
                                    if (auto f = self(newCur, child)) return true;
                                }
                                break;
                            }
                        }
                    }
                }
                if (endIdx2 > startIdx2) {
                    if (!errorMatch2) {
                        // Optimization: if only matches exist
                        // do not expand cursor to all possibilities
                        for (size_t i{startIdx2}; i < endIdx2; ++i) {
                            auto const& [child1, child2] = This->stackNode2[i];
                            auto symb = std::array<size_t, KStep>{ //!TODO this only works for KStep == 2
                                child1.side[node.dir].lastIRank,
                                child2.side[node.dir].lastIRank
                            };


                            bool insertHappened = child1.side[node.dir].queryPos + dirStep != child2.side[node.dir].queryPos;
                            if (insertHappened) {
                                auto nextSymb2 = child2.side[node.dir].lastQRank;
                                auto posI = child2.side[node.dir].queryPos - dirStep*2;
                                auto prevSymbI = This->query[posI];
                                if (nextSymb2 == prevSymbI) continue;
                            }

                            // match
                            auto newCur = This->extendKStep(cur, symb, node.dir);
                            if (newCur.count() == 0) continue;
                            if (auto f = self(newCur, child2)) return true;
                        }
                    } else {
                        // Expand cursor to all k-step paths
    //                    auto cursors = This->extend(cur, node.dir);
                        auto cursors2 = This->extendKStep(cur, node.dir); //!TODO could also produce them for 1 steps
                        for (size_t i{startIdx2}; i < endIdx2; ++i) {

     /*                       auto& [child, child2] = This->stackNode2[i];
                            switch(child.side[node.dir].info) {
                            case 'M': {
                                auto const& newCur = cursors[child.side[node.dir].lastIRank];
                                if (newCur.count() == 0) continue;
                                if (auto f = self(newCur, child)) return true;
                                break;
                            }
                            }*/

                            #if 1
                            auto& [child1, child2] = This->stackNode2[i];

                            auto step2 = [&]() -> bool {
                                switch(child2.side[node.dir].info) {
                                case 'M': {
                                    auto tinfo1 = child1.side[node.dir].info;
                                    auto nextSymb2 = child2.side[node.dir].lastQRank;
                                    auto posI = child2.side[node.dir].queryPos - dirStep*2;
                                    auto prevSymbI = This->query[posI];
                                    //bool insertHappened = child1.side[node.dir].queryPos + dirStep != child2.side[node.dir].queryPos;
                                    //Here is the mistake, tinfo1 is not the previous node, but the one before that (last 'M' movement)
                                    //This causes do let this one falsly slide through.....
                                    bool insertHappened = child1.side[node.dir].queryPos + dirStep != child2.side[node.dir].queryPos;
                                    if (insertHappened) {
                                        auto nextSymb2 = child2.side[node.dir].lastQRank;
                                        auto posI = child2.side[node.dir].queryPos - dirStep*2;
                                        auto prevSymbI = This->query[posI];
                                        if (nextSymb2 == prevSymbI) break;;
                                    }
    //                                if (tinfo1 == 'I' && nextSymb2 == child1.side[node.dir].lastQRank) break;
                                    if (tinfo1 == 'D' && nextSymb2 == child1.side[node.dir].lastIRank) break;
                                    auto symb1 = child1.side[node.dir].lastIRank;
                                    auto symb2 = child2.side[node.dir].lastIRank;
                                    auto symb = symb1 * Sigma + symb2;
                                    auto const& newCur = cursors2[symb];
                                    if (newCur.count() == 0) return false;;
                                    if (auto f = self(newCur, child2)) return true;
                                    return false;
                                }
                                case 'S': {
                                    auto symb1 = child1.side[node.dir].lastIRank;
                                    for (uint64_t i2{FirstSymb}; i2 < Sigma; ++i2) {
                                        if (i2 == child2.side[node.dir].lastQRank) continue;
                                        auto symb = symb1 * Sigma + i2;
                                        auto const& newCur = cursors2[symb];
                                        if (newCur.count() == 0) continue;
                                        child2.side[node.dir].lastIRank = i;
                                        if (auto f = self(newCur, child2)) return true;
                                    }
                                    return false;
                                }
                                case 'D': {
                                    auto symb1 = child1.side[node.dir].lastIRank;
                                    for (uint64_t i2{FirstSymb}; i2 < Sigma; ++i2) {
                                        auto symb = symb1 * Sigma + i2;
                                        auto const& newCur = cursors2[symb];
                                        if (newCur.count() == 0) continue;
                                        child2.side[node.dir].lastIRank = i2;
                                        if (auto f = self(newCur, child2)) return true;
                                    }
                                    return false;
                                }
                                case 'I': {
    //                                if (auto f = self(cur, child2)) return true;
                                    throw std::runtime_error{"this should not happen 1"};
                                    //!TODO
                                }
                                }
                                return false;
                            };

                            switch(child1.side[node.dir].info) {
                            case 'M': {
                                if (step2()) return true;
                                break;
                            }
                            case 'S':
                                for (uint64_t i{FirstSymb}; i < Sigma; ++i) {
                                    if (i == child1.side[node.dir].lastQRank) continue;
                                    child1.side[node.dir].lastIRank = i;
                                    if (step2()) return true;
                                }
                                break;
                            case 'D':
                                for (uint64_t i{FirstSymb}; i < Sigma; ++i) {
                                    child1.side[node.dir].lastIRank = i;
                                    if (step2()) return true;
                                }
                                break;
                            case 'I':
                                throw std::runtime_error{"does not happen"};
                            }
                        #endif
                        }
                    }
                }

                for (auto const& child : optInsert) {
                    // inserted into query
                    if (auto f = self(cur, child)) return true;
                }
                if (optInsert2.size() > 0) {
                    auto cursors = This->extend(cur, node.dir);
                    for (auto const& [child1, _child2] : optInsert2) {
                        auto child2 = _child2;
                        switch(child1.side[node.dir].info) {
                        case 'M': {
                            auto const& newCur = cursors[child1.side[node.dir].lastIRank];
                            if (newCur.count() == 0) continue;
                            //if (auto f = self(newCur, child1)) return true;
                            if (auto f = self(newCur, child2)) return true;
                            break;
                        }
                        case 'S':
                            for (uint64_t i{FirstSymb}; i < Sigma; ++i) {
                                if (i == child1.side[node.dir].lastQRank) continue;
                                auto const& newCur = cursors[i];
                                if (newCur.count() == 0) continue;
                                child2.side[node.dir].lastIRank = i;
                                if (auto f = self(newCur, child2)) return true;
                            }
                            break;
                        case 'D':
                            for (uint64_t i{FirstSymb}; i < Sigma; ++i) {
                                auto const& newCur = cursors[i];
                                if (newCur.count() == 0) continue;
                                child2.side[node.dir].lastIRank = i;
                                if (auto f = self(newCur, child2)) return true;
                            }
                            break;
                        case 'I':
                            throw std::runtime_error{"does not happen"};
                        }
                    }
                }
                /*if (optInsert) {
                    auto const& child = *optInsert;
                    // inserted into query
                    if (auto f = self(cur, child)) return true;
                }*/
                This->stackNode.resize(startIdx);
                This->stackNode2.resize(startIdx2);
            } else {
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

                if (!errorMatch) {
                    // Optimization: if only matches exist
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
                        switch(child.side[node.dir].info) {
                        case 'M': {
                            auto const& newCur = cursors[child.side[node.dir].lastIRank];
                            if (newCur.count() == 0) continue;
                            if (auto f = self(newCur, child)) return true;
                            }
                            break;
                        case 'S':
                            for (uint64_t i{FirstSymb}; i < Sigma; ++i) {
                                if (i == child.side[node.dir].lastQRank) continue;
                                auto const& newCur = cursors[i];
                                if (newCur.count() == 0) continue;
                                child.side[node.dir].lastIRank = i;
                                if (auto f = self(newCur, child)) return true;
                            }
                            break;
                        case 'D':
                            for (uint64_t i{FirstSymb}; i < Sigma; ++i) {
                                auto const& newCur = cursors[i];
                                if (newCur.count() == 0) continue;
                                child.side[node.dir].lastIRank = i;
                                if (auto f = self(newCur, child)) return true;
                            }
                            break;
                        }
                    }
                }

                if (optInsert) {
                    auto const& child = *optInsert;
                    // inserted into query
                    if (auto f = self(cur, child)) return true;
                }
                This->stackNode.resize(startIdx);
            }

            return false;
        };
        auto cur = cursor_t{index};
        return handleNode(cur, rootNode);
    }

    void expandNode_2Step(Node node, auto const& cb1, auto const& cb2) {
        expandNode(node, [&](Node child1) {
            // end of parts, direction change might happen or end of search
            // ?TODO this could be combined with the next one and a check if direction has changed could be performed
            if (child1.partitionPart == 0) {
                cb1(child1);
                return;
            }
/*            cb1(child1);
            return;*/
            switch (child1.side[node.dir].info) {
            case 'M':
                expandNode<true>(child1, [&](Node const& child2) {
//                    if (child2.side[node.dir].info == 'I') throw std::runtime_error{"this should not happen 2"};
/*                    if (child2.side[node.dir].info == 'I') {
                        cb1(child1);
                    } else*/ {
                        cb2(child1, child2);
                    }
                });
                break;
            #if 1
            case 'S':
                expandNode<true>(child1, [&](Node const& child2) {
/*                    if (child2.side[node.dir].info == 'I') {
                        cb1(child1);
                    } else*/ {
                        cb2(child1, child2);
                    }
                });
                break;
            case 'D':
                expandNode<true>(child1, [&](Node const& child2) {
/*                    if (child2.side[node.dir].info == 'I') {
                        cb1(child1);
                    } else*/ {
                        cb2(child1, child2);
                    }
                });
                break;
            #endif
            case 'I': // should be covered by the previous partitionPart == 0 case
                throw std::runtime_error{"this should not happen 3"};
                cb1(child1); // This kind of implies it
                break;
            }
        });
    }
};


template <bool Edit, Optimizations OPT, typename index_t, Sequence query_t, typename report_t>
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
        bool f = Search<Edit, index_t, query_t, decltype(search), decltype(internal_report), OPT>{index, query, search, partition, internal_report}.run();
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
template <bool Edit, Optimizations OPT, typename index_t, Sequences queries_t, typename selectSearchScheme_t, typename report_t>
void search_n_impl(index_t const& index, queries_t&& queries, selectSearchScheme_t&& selectSearchScheme, std::vector<size_t> const& partition, report_t&& report, size_t n) {
    if (queries.empty()) return;
    if (n == 0) return;
    for (size_t qidx{}; qidx < queries.size(); ++qidx) {
        size_t ct{};
        auto const& search_scheme = selectSearchScheme(queries[qidx].size());
        search_impl<Edit, OPT>(index, queries[qidx], search_scheme, partition, [&] (auto cur, size_t e) {
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
template <bool Edit=true, Optimizations OPT = {}, typename index_t, Sequences queries_t, typename report_t>
void search(index_t const& index, queries_t&& queries, search_scheme::Scheme const& search_scheme, std::vector<size_t> const& partition, report_t&& report, size_t n = std::numeric_limits<size_t>::max()) {
    // function that selects a search scheme
    auto selectSearchScheme = [&]([[maybe_unused]] size_t length) -> auto& {
        return search_scheme;
    };
    search_n_impl<Edit, OPT>(index, queries, selectSearchScheme, partition, report, n);
}

}
