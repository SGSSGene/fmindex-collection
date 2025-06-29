// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <cstdint>
#include <optional>
#include <tuple>
#include <vector>

namespace fmc {


template <typename index_t, typename cursor_t>
struct LocateLinear {
    struct iter {
        index_t const& index;
        size_t pos;

        auto operator*() const -> std::tuple<typename index_t::ADEntry, size_t> {
            return index.locate(pos);
        }
        auto operator++() -> iter& {
            pos += 1;
            return *this;
        }
        auto operator!=(size_t _end) const {
            return pos != _end;
        }
    };
    //! same as iter, but for reversed text fm index
    struct rev_iter : iter {
        size_t depth;

        auto operator*() const -> std::tuple<typename index_t::ADEntry, size_t> {
            auto [entry, offset] = iter::index.locate(iter::pos);
            assert(offset >= depth);
            return {entry, offset - depth};
        }
    };


    index_t const& index;
    cursor_t cursor;

    friend auto begin(LocateLinear const& locate) {
        // Check if it is the reversed cursor
        if constexpr (requires(cursor_t c) {{c.query_length()} -> std::same_as<size_t>; }) {
            return rev_iter{locate.index, locate.cursor.lb, locate.cursor.depth};
        } else {
            return iter{locate.index, locate.cursor.lb};
        }
    }
    friend auto end(LocateLinear const& locate) {
        return locate.cursor.lb + locate.cursor.len;
    }
};

//!TODO remove as soon as clang supports auto deduction guides (not the case in clang 15
#if __clang__
    template <typename index_t, typename cursor_t>
    LocateLinear(index_t const&, cursor_t) -> LocateLinear<index_t, cursor_t>;
#endif

template <typename index_t, typename cursor_t>
struct LocateFMTree {
    std::vector<std::tuple<typename index_t::ADEntry, size_t>> positions;

    LocateFMTree(index_t const& index, cursor_t cursor, size_t samplingRate, size_t maxDepth) {
        positions.reserve(cursor.count());
        auto stack = std::vector<std::tuple<cursor_t, size_t>>{};
        stack.emplace_back(cursor, 0);
        while(!stack.empty()) {
            auto [cursor, depth] = stack.back();
            stack.pop_back();
            //!TODO what should the termination criteria be?
            if (depth >= maxDepth or cursor.count() < index_t::Sigma*2) {
                for (size_t pos{cursor.lb}; pos < cursor.lb + cursor.len; ++pos) {
                    uint64_t maxSteps = samplingRate - depth;
                    uint64_t idx = pos;

                    auto opt = index.single_locate_step(idx);
                    uint64_t steps{};
                    for (;!opt && maxSteps > 0; --maxSteps) {
                        auto symb = index.bwt.symbol(idx);
                        idx = index.bwt.rank(idx, symb) + index.C[symb];
                        steps += 1;
                        opt = index.single_locate_step(idx);
                    }
                    if (opt) {
                        positions.emplace_back(*opt, depth+steps);

                    }
                }
            } else {
                for (size_t pos{cursor.lb}; pos < cursor.lb + cursor.len; ++pos) {
                    auto v = index.single_locate_step(pos);
                    if (v) {
                        positions.emplace_back(*v, depth);
                    }
                }

                auto cursors = cursor.extendLeft();
                for (size_t sym{1}; sym < cursors.size(); ++sym) {
                    stack.emplace_back(cursors[sym], depth+1);
                }
            }
        }
        assert(positions.size() == cursor.count());
    }

    friend auto begin(LocateFMTree const& locate) {
        return begin(locate.positions);
    }
    friend auto end(LocateFMTree const& locate) {
        return end(locate.positions);
    }
};


template <size_t MaxDepth, typename index_t, typename cursor_t, typename CB>
void locateFMTree(index_t const& index, cursor_t cursor, CB const& cb, size_t samplingRate, size_t depth=0) {
    if (depth < MaxDepth and cursor.count() > 1000) {
        for (size_t pos{cursor.lb}; pos < cursor.lb + cursor.len; ++pos) {
            auto v = index.single_locate_step(pos);
            if (v) {
                cb(*v, depth);
            }
        }

        if (depth+1 < samplingRate) {
            auto cursors = cursor.extendLeft();
            for (size_t sym{1}; sym < cursors.size(); ++sym) {
                locateFMTree<MaxDepth>(index, cursors[sym], cb, depth+1);
            }
        }
    } else {
        for (size_t pos{cursor.lb}; pos < cursor.lb + cursor.len; ++pos) {
            uint64_t maxSteps = samplingRate - depth;
            uint64_t idx = pos;

            auto opt = index.single_locate_step(idx);
            uint64_t steps{};
            for (;!opt && maxSteps > 0; --maxSteps) {
                auto symb = index.bwt.symbol(idx);
                idx = index.bwt.rank(idx, symb) + index.C[symb];
                steps += 1;
                opt = index.single_locate_step(idx);
            }
            if (opt) {
                cb(*opt, depth + steps);
            }
        }
    }
}

}
