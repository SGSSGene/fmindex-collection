#pragma once

#include <vector>
#include <cstddef>
#include <cstdint>

#include <fmt/format.h>


/**
 * Ng16 with replaced bandedmatrix
 */
namespace search_ng17 {

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
        data.resize(rows*cols);
    }

    size_t operator[](size_t idx) const {
        return data[idx];
    }

    size_t& operator[](size_t idx) {
        return data[idx];
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
        matrix[0] = e;
        search_next(_cursor, 0, 0, 1);
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

    void search_next(cursor_t const& cur, size_t pos, size_t start, size_t end) noexcept {
        auto cursors = extend(cur);

        for (size_t symb{1}; symb < Sigma; ++symb) {
            auto const& cur = cursors[symb];
            if (cur.empty()) continue;

            auto length = end-start;

            auto newPos   = pos;
            auto newStart = end;
            auto newEnd   = end+length;
            // init first cell
            {
                auto const& s = search[pos];
                auto rank = query[s.pi];

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
            // fill last match/mismatch
            if (pos + length < search.size()) {
                newEnd += 1;
                auto const& s = search[pos+length];
                auto rank = query[s.pi];

                auto gapX = matrix[newStart+length-1]+1;
                auto diag = matrix[start+length-1] + (rank != symb);

                matrix[newStart+length] = std::min({gapX, diag});
            }
            // fill extra
            {
                auto isValid = [&](size_t pos, size_t i) {
                    if (pos >= search.size()) return false;
                    auto const& s = search[pos];
                    auto rank = query[s.pi];

                    auto e = matrix[i];
                    return s.l <= e and e <= s.u;
                };
                bool valid = isValid(length, newEnd-1);
                size_t i{1};
                while (valid) {
                    matrix[newEnd] = matrix[newEnd-1]+1;
                    valid = isValid(length+i, newEnd);
                    newEnd += 1;
                    i += 1;
                }
            }

            // find valid start
            while (newStart < newEnd) {
                auto const& s = search[newPos];

                auto e = matrix[newStart];
                if (s.l <= e and e <= s.u) {
                    break;
                }
                newPos += 1;
                newStart += 1;
            };
            // find last valid start
            for (size_t i{0}; i < newEnd - newStart; ++i) {
                auto const& s = search[newPos+i];

                auto e = matrix[newStart+i];
                if (!(s.l <= e and e <= s.u)) {
                    newEnd = newStart+i;
                    break;
                }
            };

            if (newStart < newEnd) {
                if (newPos + (newEnd-newStart) == search.size()) {
                    delegate(cur, matrix[newEnd-1]);
                }
                if (newPos < search.size()) {
                    search_next(cur, newPos, newStart, newEnd);
                }
            }
        }
    }
};



template <typename index_t, typename queries_t, typename search_schemes_t, typename delegate_t>
void search(index_t const & index, queries_t && queries, search_schemes_t const & search_scheme, delegate_t && delegate)
{
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
            search2.back().emplace_back(Block{(size_t)s.pi[i], (size_t)s.l[i], (size_t)s.u[i]});
        }
        search_scheme2.emplace_back(move(search2));
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
