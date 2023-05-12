// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "../BiFMIndexCursor.h"

#include <vector>
#include <cstddef>
#include <cstdint>

/*
 * using banded alignment matrix, but only works for backtracking search schemes
 */
namespace fmindex_collection {
namespace search_ng20 {

struct BandMatrix {
private:
    size_t queryLength;
    size_t maxErrors;

    size_t rows;
    size_t cols;
    std::vector<size_t> data;
public:

    BandMatrix(size_t _query, size_t _maxErrors)
        : queryLength{_query}
        , maxErrors{_maxErrors}
        , rows{queryLength+maxErrors+1}
        , cols{queryLength+1}
    {
        data.resize(rows*cols, maxErrors + 1);
        for (size_t i{0}; i <= std::min(cols, _maxErrors+1); ++i) {
            at(0, i) = i;
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
            if (row < maxErrors) {
                return size_t{0};
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

template <typename index_t, typename search_scheme_t, typename delegate_t>
struct Search {
    constexpr static size_t Sigma = index_t::Sigma;

    using cursor_t = BiFMIndexCursor<index_t>;

    index_t const& index;
    decltype(search_scheme_t::pi) const& pi;
    decltype(search_scheme_t::l) const& l;
    decltype(search_scheme_t::u) const& u;

    using query_t = std::vector<uint8_t>;


    query_t const& query;

    delegate_t const& delegate;
    size_t maxError;

    BandMatrix& matrix;

    Search(index_t const& _index, search_scheme_t const& _search, query_t const& _query, delegate_t const& _delegate, size_t _maxError, BandMatrix& _matrix) noexcept
        : index {_index}
        , pi{_search.pi}
        , l{_search.l}
        , u{_search.u}
        , query{_query}
        , delegate  {_delegate}
        , maxError{_maxError}
        , matrix{_matrix}
//        , matrix(_search.pi.size(), maxError)
    {
        auto cur = cursor_t{index};

//        search_first_allowed_error<true>(cur, 0);
        search<true>(cur, 1);
    }

    auto extend(cursor_t const& cur, uint8_t symb, std::size_t pos) const noexcept {
        if (pos == 0 or pi[pos-1] < pi[pos]) {
            return cur.extendRight(symb);
        } else {
            return cur.extendLeft(symb);
        }
    }
    template <bool Right>
    auto extend(cursor_t const& cur, uint8_t symb) const noexcept {
        if constexpr(Right) {
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

    template <bool Right>
    void search_first_allowed_error(cursor_t const& cur, size_t pos) noexcept {
        if (cur.empty()) return;

        if (pos == query.size()) {
            delegate(cur, 0);
            return;
        }

        if (u[pos] > 0) {
            for (size_t i{0}; i <= maxError+1; ++i) {
                matrix.at(pos, pos+i) = i;
            }

            search<Right>(cur, pos+1);
            return;
        }

        auto rank = query[pi[pos]];
        auto newCur = extend<Right>(cur, rank);
        search_first_allowed_error<Right>(newCur, pos+1);
    }

    template <bool Right>
    void search_error_free(cursor_t const& cur, size_t pos) noexcept {
        if (cur.empty()) return;

        if (pos == query.size()) {
            delegate(cur, 0);
            return;
        }

        auto rank = query[pi[pos]];
        auto newCur = extend<Right>(cur, rank);
        search_error_free<Right>(newCur, pos+1);
    }


    template <bool Right>
    void search(cursor_t const& cur, size_t depth) noexcept {


        auto [start, end] = matrix.rangeCol(depth);
        auto cursors = extend<Right>(cur);

        for (size_t symb{1}; symb < Sigma; ++symb) {
            auto const& cur = cursors[symb];
            if (cur.empty()) continue;

            bool atLeastOneInside = false;
            size_t smallestValidError = std::numeric_limits<size_t>::max();
            for (size_t i{start}; i < end; ++i) {
                auto rank = query[pi[i-1]];
                auto e = matrix.updateCell(depth, i, rank == symb);
                if (l[i-1] <= e and e <= u[i-1]) {
                    smallestValidError = std::min(smallestValidError, e);
                    atLeastOneInside = true;
                }
            }
            if (end == query.size() + 1) {
                auto e = matrix.at(depth, end-1);
                if (l[end-2] <= e and e <= u[end-2]) {
                    delegate(cur, e);
                }
            }
            if (atLeastOneInside) {
                if (smallestValidError == maxError) {
                    for (size_t i{start}; i < end; ++i) {
                        auto e = matrix.at(depth, i);
                        if (l[i-1] <= e and e <= u[i-1]) {
                            search_error_free<Right>(cur, i);
                        }
                    }
/*                if (end == query.size()+1) {
                    delegate(cur, smallestValidError);
//                    fmt::print("delegate: {}\n", depth);
//                    matrix.print(depth);
//                    fmt::print("\n");
                    continue;
                }*/
                } else {
                    search<Right>(cur, depth+1);
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

    std::size_t qidx;
    auto internal_delegate = [&qidx, &delegate] (auto const & it, size_t e) {
        delegate(qidx, it, e);
    };
    size_t maxError{0};
    for (auto const& s : search_scheme) {
        for (auto const& u : s.u) {
            maxError = std::max(maxError, size_t(u));
        }
    }

    BandMatrix matrix{queries[0].size(), maxError};

//    for (qidx = {0}; qidx < 1; ++qidx) {
    for (qidx = {0}; qidx < queries.size(); ++qidx) {

        for (size_t j{0}; j < search_scheme.size(); ++j) {
            Search<std::decay_t<decltype(index)>, std::decay_t<decltype(search_scheme[j])>, std::decay_t<decltype(internal_delegate)>> {index, search_scheme[j], queries[qidx], internal_delegate, maxError, matrix};
        }
    }
}

}
}
