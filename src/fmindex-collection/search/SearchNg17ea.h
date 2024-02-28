// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "SelectCursor.h"

#include <vector>
#include <cstddef>
#include <cstdint>
#include <functional>

/**
 * Like Ng17, but aborts early
 * !TODO intext verification, configurable early abort
 */
namespace fmindex_collection::search_ng17ea {

struct BandMatrix {
private:
public:
    size_t queryLength;
    size_t maxErrors;

    size_t rows;
    size_t cols;
    std::vector<size_t> data;

    struct View {
        size_t* data;
        size_t size;

        size_t operator[](size_t idx) const {
            assert(idx < size);
            return data[idx];
        }
        size_t& operator[](size_t idx) {
            assert(idx < size);
            return data[idx];
        }
        void push(size_t v) {
            data[size] = v;
            ++size;
        }
        size_t back() const {
            assert(size > 0);
            return data[size-1];
        }
        size_t front() const {
            assert(size > 0);
            return data[0];
        }
        void pop_front() {
            assert(size > 0);
            data = data+1;
            --size;
        }
    };
public:

    BandMatrix(size_t _query, size_t _maxErrors)
        : queryLength{_query}
        , maxErrors{_maxErrors}
        , rows{queryLength+maxErrors+1}
        , cols{queryLength+1}
    {
        data.resize(rows*cols);
    }

    size_t operator[](size_t idx) const {
        return data[idx];
    }

    size_t& operator[](size_t idx) {
        return data[idx];
    }

    auto row(size_t idx, size_t size) {
        return View{data.data()+idx, size};
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
        , matrix{_search.size(), _maxError}
    {
        matrix[0] = e;
        size_t newEnd = 1;
        while(newEnd < search.size() && matrix[newEnd-1]+1 <= search[newEnd].u) {
            matrix[newEnd] = matrix[newEnd-1] + 1;
            ++newEnd;
        }

        search_next2(_cursor, 0, 0, newEnd);
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

        if (cur.count() < search.size() - pos) {
            delegate(cur, e);
            return;
        }

        if (pos == search.size()) {
            delegate(cur, e);
            return;
        }

        auto rank = query[search[pos].pi];
        auto newCur = extend(cur, rank);
        search_error_free(newCur, e, pos+1);
    }


    void search_next2(cursor_t const& cur, size_t const pos, size_t const start, size_t end) noexcept {
        auto length = end-start;

        if (length == 0) {
            return;
        }
        if (cur.empty()) {
            return;
        }
        if (cur.count() < search.size() - pos) {
            delegate(cur, matrix[end-1]);
            return;
        }


        if (pos + length == search.size()) {
            delegate(cur, matrix[end-1]);
            end -= 1;
            length -= 1;
            if (length == 0) {
                return;
            }
        }

        auto cursors = extend(cur);

        assert(pos + length <= search.size());

        for (size_t symb{1}; symb < Sigma; ++symb) {
            auto const& cur = cursors[symb];
            if (cur.empty()) continue;


            auto newPos   = pos;
            auto newStart = end;
            auto newEnd   = newStart+length;
            // init first cell
            {
                matrix[newStart] = matrix[start]+1;
            }

            // fill others
            for (size_t i{1}; i < length; ++i) {
                auto const& s = search[pos+i];
                auto rank = query[s.pi];

                auto gapX = matrix[newStart+i-1]+1;
                auto gapY = matrix[start+i]+1;
                auto diag = matrix[start+i-1] + (rank != symb);

                matrix[newStart+i] = std::min({gapX, gapY, diag});
            }

            if (newPos + length < search.size()) {
                // fill last match/mismatch
                {
                    auto const& s = search[newPos + length];
                    auto rank = query[s.pi];

                    auto gapX = matrix[newEnd-1]+1;
                    auto diag = matrix[start+length-1] + (rank != symb);

                    matrix[newEnd] = std::min({gapX, diag});
                    newEnd += 1;
                }
            }

            // find valid start
            {
                while (newStart < newEnd) {
                    assert(newPos < search.size());
                    auto const& s = search[newPos];

                    auto e = matrix[newStart];
                    if (s.l <= e and e <= s.u) {
                        break;
                    }
                    newPos += 1;
                    newStart += 1;
                };
            }
            assert(newPos + newEnd - newStart <= search.size());
            // fill extra
            {
                auto isValid = [&](size_t offset, size_t i) {
                    if (newPos + offset >= search.size()) return false;
                    auto const& s = search[newPos + offset];

                    auto e = matrix[i]+1;
                    return s.l <= e and e <= s.u;
                };

                size_t i{1};
                bool valid = isValid(length+i, newEnd-1);
                while (valid) {
                    matrix[newEnd] = matrix[newEnd-1]+1;
                    newEnd += 1;
                    i += 1;
                    valid = isValid(length+i, newEnd-1);
                }
            }

            // invalidate, invalid entries
            auto newLength = newEnd - newStart;
            for (size_t i{0}; i < newLength; ++i) {
                assert(newPos + newLength - i - 1 < search.size());

                auto const& s = search[newPos+newLength-i-1];
                auto& e = matrix[newStart+newLength-i-1];
                if ((s.l <= e and e <= s.u)) {
                    newEnd = newEnd - i;
                    break;
                }

            }
            /*for (size_t i{0}; i < newEnd - newStart; ++i) {
                assert(newPos+i < search.size());
                auto const& s = search[newPos+i];

                auto& e = matrix[newStart+i];
                if (!(s.l <= e and e <= s.u)) {
                    e = std::numeric_limits<size_t>::max()/2;
//                    newEnd = newStart+i;
//                    break;
                }
            };*/


            if (newStart < newEnd) {
                search_next2(cur, newPos, newStart, newEnd);
            }
        }
    }

    void search_next(cursor_t const& cur, size_t const pos, size_t const start, size_t end) noexcept {
        if (cur.empty()) {
            return;
        }
        if (cur.count() < search.size() - pos) {
            delegate(cur, matrix[end-1]);
            return;
        }

        auto length = end-start;

        if (pos + length == search.size()) {
            delegate(cur, matrix[end-1]);
        }
        end = std::min(end, start + (search.size() - pos));
        length = end - start;
        if (pos >= search.size() or length == 0) {
            return;
        }
        auto cursors = extend(cur);

        assert(pos + length <= search.size());

        auto const lastRow = matrix.row(start, end-start);

        for (size_t symb{1}; symb < Sigma; ++symb) {
            auto const& cur = cursors[symb];
            if (cur.empty()) continue;

            auto newPos = pos;
            auto newRow = matrix.row(end, 0);

            // init first cell
            {
                newRow.push(lastRow[0]+1);
            }
            // fill others
            for (size_t i{1}; i < length; ++i) {
                auto const& s = search[pos+i];
                auto rank = query[s.pi];

                auto gapX = newRow.back()+1;
                auto gapY = lastRow[i]+1;
                auto diag = lastRow[i-1] + (rank != symb);

                newRow.push(std::min({gapX, gapY, diag}));
            }

            if (newPos + newRow.size < search.size()) {
                // fill last match/mismatch
                {
                    auto const& s = search[newPos + newRow.size];
                    auto rank = query[s.pi];

                    auto gapX = newRow.back()+1;
                    auto diag = lastRow.back() + (rank != symb);

                    newRow.push(std::min({gapX, diag}));
                }

                // fill extra
                {
                    auto isValid = [&]() {
                        if (newPos + newRow.size >= search.size()) return false;
                        auto const& s = search[newPos + newRow.size];

                        auto e = lastRow.back();
                        return s.l <= e and e <= s.u;
                    };

                    while (isValid()) {
                        newRow.push(lastRow.back()+1);
                    }
                }
            }
            // find valid start
            {
                while (newRow.size > 0) {
                    assert(newPos < search.size());
                    auto const& s = search[newPos];

                    auto e = newRow.front();
                    if (s.l <= e and e <= s.u) {
                        break;
                    }
                    newPos += 1;
                    newRow.pop_front();
                };
            }
            // find last valid start
            for (size_t i{0}; i < newRow.size; ++i) {
                assert(newPos+i < search.size());
                auto const& s = search[newPos+i];

                auto e = newRow[i];
                if (!(s.l <= e and e <= s.u)) {
                    newRow.size = i;
                    break;
                }
            };

            if (newRow.size > 0) {
                auto newStart = newRow.data - matrix.data.data();
                auto newEnd   = newStart + newRow.size;
                search_next(cur, newPos, newStart, newEnd);
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
                search2.back().emplace_back(Block{std::numeric_limits<size_t>::max(), size_t{s.l[i]}, size_t{s.u[i]}});
            } else if (lastDir != dir) {
                auto lastEntry = search2.back().back();
                lastEntry.pi = std::numeric_limits<size_t>::max();
                search2.emplace_back();
                search2.back().emplace_back(lastEntry);

            }
            if (search2.back().back().pi == std::numeric_limits<size_t>::max()) {
                search2.back().back().u = size_t{s.u[i]};
            }
            search2.back().emplace_back(Block{size_t{s.pi[i]}, size_t{s.l[i]}, size_t{s.u[i]}});
            lastDir = dir;

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
        if (search->at(depth).empty()) {
            depth += 1;
            sch(cursor, e);
            depth -= 1;
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
