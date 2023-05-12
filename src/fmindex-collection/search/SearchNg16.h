// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "../BiFMIndexCursor.h"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <vector>

/**
 * As ng15 and 20 combined
 */
namespace fmindex_collection {
namespace search_ng16 {

struct BandMatrix {
private:
    size_t queryLength;
    size_t maxErrors;

    size_t rows;
    size_t cols;
    std::vector<size_t> data;
public:

    BandMatrix(size_t _query, size_t e, size_t _maxErrors)
        : queryLength{_query}
        , maxErrors{_maxErrors}
        , rows{queryLength+maxErrors+1}
        , cols{queryLength+1}
    {
        data.resize(rows*cols, maxErrors + 1);
        for (size_t i{0}; i <= std::min(cols, _maxErrors+1); ++i) {
            at(0, i) = e+i;
        }
    }

    size_t at(size_t row, size_t col) const {
        auto [begin, end] = rangeCol(row);
        return data[row * (maxErrors+3) + col - begin+1];
    }
    size_t& at(size_t row, size_t col) {
        auto [begin, end] = rangeCol(row);
        return data[row * (maxErrors+3) + col - begin+1];
//        return data[row * cols + col];
    }

    size_t updateCell(size_t row, size_t col, bool match) {
        auto diag = at(row-1, col-1) + !match; // match/mismatch
        auto gapY = at(row-1, col) + 1;        // substitution
        auto gapX = at(row, col - 1) + 1;      // deletion
        auto v = std::min({diag, gapY, gapX});
        at(row, col) = v;
        return v;
    }

    auto rangeCol(size_t row) const -> std::tuple<size_t, size_t> {
        auto start = [&] () {
//            return row;
            if (row <= maxErrors) {
                return size_t{1};
            } else {
                return row - maxErrors;
            }
        }();
        auto end = [&] () {
//            return std::min(cols, row+1);
            if (row + maxErrors+1 > cols) {
                return cols;
            } else {
                return row + maxErrors+1;
            }
        }();
        return {start, end};
    }
};


struct Block {
    size_t pi;
    size_t l;
    size_t u;
};

template <typename cursor_t, typename search_scheme_t, typename query_t, typename delegate_t, bool Right>
struct Search {
    constexpr static size_t Sigma = cursor_t::Sigma;

    search_scheme_t const& search;
    query_t const& query;
    delegate_t const& delegate;

    size_t maxError;

    BandMatrix matrix;


    Search(cursor_t const& _cursor, search_scheme_t const& _search, query_t const& _query, size_t e, delegate_t const& _delegate, size_t _maxError)
        : search    {_search}
        , query     {_query}
        , delegate  {_delegate}
        , maxError{_maxError}
        , matrix{_search.size(), e, _maxError}
    {
        search_next(_cursor, 1);
    }

    static auto extend(cursor_t const& cur, uint8_t symb) noexcept {
        if constexpr (Right) {
            return cur.extendRight(symb);
        } else {
            return cur.extendLeft(symb);
        }
    }
    static auto extend(cursor_t const& cur) noexcept {
        if constexpr (Right) {
            return cur.extendRight();
        } else {
            return cur.extendLeft();
        }
    }

    void search_error_free(cursor_t const& cur, size_t e, size_t pos) noexcept {
        if (cur.empty()) return;

        if (pos == search.size()) {
            delegate(cur, e);
            return;
        }

        auto rank = query[search[pos].pi];
        auto newCur = extend(cur, rank);
        search_error_free(newCur, e, pos+1);
    }

    void search_next(cursor_t const& cur, size_t depth) noexcept {
        auto [start, end] = matrix.rangeCol(depth);
        auto cursors = extend(cur);

        for (size_t symb{1}; symb < Sigma; ++symb) {
            auto const& cur = cursors[symb];
            if (cur.empty()) continue;

            bool atLeastOneInside = false;
            size_t smallestValidError = std::numeric_limits<size_t>::max();
            for (size_t i{start}; i < end; ++i) {
                auto rank = query[search[i-1].pi];
                auto e = matrix.updateCell(depth, i, rank == symb);
                if (search[i-1].l <= e and e <= search[i-1].u) {
                    smallestValidError = std::min(smallestValidError, e);
                    atLeastOneInside = true;
                }
            }
            if (end == search.size() + 1) {
                auto e = matrix.at(depth, end-1);
                if (search[end-2].l <= e and e <= search[end-2].u) {
                    delegate(cur, e);
                }
            }
            if (atLeastOneInside) {
                if (smallestValidError == maxError) {
                    for (size_t i{start}; i < end; ++i) {
                        auto e = matrix.at(depth, i);
                        if (search[i-1].l <= e and e <= search[i-1].u) {
                            search_error_free(cur, e, i);
                        }
                    }
                } else {
                    search_next(cur, depth+1);
                }
            } else {
//                fmt::print("out of range: {}\n", depth);
//                matrix.print(depth);
//                fmt::print("\n");

            }
        }
    }
};



template <typename index_t, typename queries_t, typename search_schemes_t, typename delegate_t>
void search(index_t const & index, queries_t && queries, search_schemes_t const & search_scheme, delegate_t && delegate)
{
    using cursor_t = select_cursor_t<index_t>;
    static_assert(not cursor_t::Reversed, "reversed fmindex is not supported");

    if (search_scheme.empty()) return;

    std::vector<std::vector<std::vector<Block>>> search_scheme2;
    size_t maxError{0};
    for (auto s : search_scheme) {
        for (auto const& u : s.u) {
            maxError = std::max(maxError, size_t(u));
        }

        std::vector<std::vector<Block>> search2;
        int lastDir = 0;
        for (size_t i{0}; i < s.pi.size(); ++i) {

            auto dir = [&]() {
                if (i == 0) {
                    return s.pi[i] < s.pi[i+1]?1:-1;
                } else {
                    return s.pi[i-1] < s.pi[i]?1:-1;
                }
            }();
            if (lastDir == 0) {
                search2.emplace_back();
                if (dir == -1) {
                    search2.emplace_back();
                }
            } else if (lastDir != dir) {
                search2.emplace_back();
            }
            lastDir = dir;
            search2.back().emplace_back(Block{size_t{s.pi[i]}, size_t{s.l[i]}, size_t{s.u[i]}});
        }
        search_scheme2.emplace_back(std::move(search2));
    }

    using Cursor = BiFMIndexCursor<index_t>;


    using Callback = std::function<void(Cursor const&, size_t e)>;
    Callback sch;

    using query_t = std::decay_t<decltype(queries[0])>;
    using search_t = std::decay_t<decltype(search_scheme2[0][0])>;

    size_t qidx;
    query_t const* query;
    std::vector<search_t> const* search;
    size_t depth{};
    sch = [&](Cursor const& cursor, size_t e) {
        if (depth == search->size()) {
            delegate(qidx, cursor, e);
            return;
        }
        if (depth % 2 == 0) {
            Search<Cursor, search_t, query_t, Callback, true>{cursor, search->at(depth++), *query, e, sch, maxError};
        } else {
            Search<Cursor, search_t, query_t, Callback, false>{cursor, search->at(depth++), *query, e, sch, maxError};
        }
        depth -= 1;
    };


    for (size_t i{0}; i < queries.size(); ++i) {
        qidx = i;
        query = &queries[qidx];
        for (size_t j{0}; j < search_scheme.size(); ++j) {
            search = &search_scheme2[j];
            // call sch
            sch(Cursor{index}, 0);
        }
    }

}

}
}
