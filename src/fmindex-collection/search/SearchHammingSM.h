// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../concepts.h"
#include "../fmindex/BiFMIndexCursor.h"
#include "../search_scheme/Scheme.h"
#include "SelectCursor.h"

/* Hamming distance search, without expanded search scheme
 * , using a scoring matrix
 */
namespace fmindex_collection::search_hamming_sm {

template <size_t QuerySigma, size_t RefSigma = QuerySigma>
struct ScoringMatrix {
    std::array<std::array<uint8_t, RefSigma>, QuerySigma> matrix;
    std::array<std::vector<size_t>, QuerySigma> noCostList;
    std::array<std::vector<size_t>, QuerySigma> costList;

    ScoringMatrix() {
        for (size_t y{1}; y < RefSigma; ++y) {
            for (size_t x{1}; x < QuerySigma; ++x) {
                if (x == y) {
                    setCost(x, y, 0);
                } else {
                    setCost(x, y, 1);
                }
            }
        }
    }

    void setCost(size_t queryRank, size_t refRank, size_t cost) {
        std::erase(noCostList[queryRank], refRank);
        std::erase(costList[queryRank], refRank);

        matrix[queryRank][refRank] = cost;
        if (cost == 1) {
            costList[queryRank].push_back(refRank);
        } else {
            noCostList[queryRank].push_back(refRank);
        }
    }
};

template <typename index_t, typename search_scheme_t, Sequence query_t, typename delegate_t, typename SM>
struct Search {
    constexpr static size_t Sigma = index_t::Sigma;

    using cursor_t = select_cursor_t<index_t>;

    index_t const& index;

    decltype(search_scheme_t::pi) const& pi;
    decltype(search_scheme_t::l) const& l;
    decltype(search_scheme_t::u) const& u;

    query_t const& query;
    SM const& sm;
    delegate_t const& delegate;
    size_t const partLen;
    size_t const extraPartLen;

    Search(index_t const& _index, search_scheme_t const& _search, query_t const& _query, SM const& _sm, delegate_t const& _delegate) noexcept
        : index {_index}
        , pi{_search.pi}
        , l{_search.l}
        , u{_search.u}
        , query{_query}
        , sm{_sm}
        , delegate  {_delegate}
        , partLen {query.size() / pi.size()}
        , extraPartLen {query.size() % pi.size()}
    {
        auto cur       = cursor_t{index};
        searchPart(cur, 0, 0);
    }

    void searchPart(cursor_t const& cur, size_t e, std::size_t part) const noexcept {
        if (cur.count() == 0) {
            return;
        }
        if (part == pi.size()) {
            if (l[part-1] <= e and e <= u[part-1]) {
                delegate(cur, e);
            }
            return;
        }
        if (e > u[part]) {
            return;
        }

        size_t start = pi[part]*partLen + std::min(extraPartLen, pi[part]);
        size_t end   = start + partLen + (pi[part] < extraPartLen ? 1 : 0);
        if (part == 0 or pi[part-1] < pi[part]) {
            searchPartRightDir(cur, e, part, start, end);
        } else {
            searchPartLeftDir(cur, e, part, start, end);
        }
    }

    void searchPartRightDir(cursor_t const& cur, size_t e, std::size_t part, size_t start, size_t end) const noexcept {
        if (cur.empty()) return;
        if (start == end) { // done with part
            if (l[part] <= e) {
                searchPart(cur, e, part+1);
            }
            return;
        }

        if (e+1 > u[part]) { // search only without errors (matches)
            searchPartRightDirNoErrors(cur, e, part, start, end);
            return;
        }

        auto rank = query[start];
        auto cursors = cur.extendRight();

        // search matches
        for (auto i : sm.noCostList[rank]) {
            searchPartRightDir(cursors[i], e, part, start+1, end);
        }
        // search substitutes
        if (e+1 <= u[part]) {
            for (auto i : sm.costList[rank]) {
                searchPartRightDir(cursors[i], e+1, part, start+1, end);
            }
        }
    }

    void searchPartRightDirNoErrors(cursor_t const& cur, size_t e, std::size_t part, size_t start, size_t end) const noexcept {
        if (cur.empty()) return;
        if (start == end) { // done with part
            searchPart(cur, e, part+1);
            return;
        }
        auto rank = query[start];
        for (auto i : sm.noCostList[rank]) {
            searchPartRightDirNoErrors(cur.extendRight(i), e, part, start+1, end);
        }
    }

    void searchPartLeftDir(cursor_t const& cur, size_t e, std::size_t part, size_t start, size_t end) const noexcept {
        if (cur.empty()) return;
        if (start == end) { // done with part
            if (l[part] <= e) {
                searchPart(cur, e, part+1);
            }
            return;
        }

        if (e+1 > u[part]) { // search only without errors (matches)
            searchPartLeftDirNoErrors(cur, e, part, start, end);
            return;
        }

        auto rank = query[end-1];
        auto cursors = cur.extendLeft();

        // search matches
        for (auto i : sm.noCostList[rank]) {
            searchPartLeftDir(cursors[i], e, part, start, end-1);
        }
        // search substitutes
        if (e+1 <= u[part]) {
            for (auto i : sm.costList[rank]) {
                searchPartLeftDir(cursors[i], e+1, part, start, end-1);
            }
        }
    }

    void searchPartLeftDirNoErrors(cursor_t const& cur, size_t e, std::size_t part, size_t start, size_t end) const noexcept {
        if (cur.empty()) return;
        if (start == end) { // done with part
            searchPart(cur, e, part+1);
            return;
        }
        auto rank = query[end-1];
        for (auto i : sm.noCostList[rank]) {
            searchPartLeftDirNoErrors(cur.extendLeft(i), e, part, start, end-1);
        }
    }


};

template <typename index_t, Sequence query_t, typename search_schemes_t, typename delegate_t, typename SM>
void search(index_t const & index, query_t && query, search_schemes_t const & search_scheme, SM const& sm, delegate_t && delegate) {
    using cursor_t = select_cursor_t<index_t>;
    static_assert(not cursor_t::Reversed, "reversed fmindex is not supported");

    for (size_t j{0}; j < search_scheme.size(); ++j) {
        Search{index, search_scheme[j], query, sm, delegate};
    }
}


template <typename index_t, Sequences queries_t, typename search_schemes_t, typename delegate_t, typename SM>
void search(index_t const & index, queries_t && queries, search_schemes_t const & search_scheme, SM const& sm, delegate_t && delegate) {
    using cursor_t = select_cursor_t<index_t>;
    static_assert(not cursor_t::Reversed, "reversed fmindex is not supported");

    std::size_t qidx;
    auto internal_delegate = [&qidx, &delegate] (auto const & it, size_t e) {
        delegate(qidx, it, e);
    };

    for (qidx = {0}; qidx < queries.size(); ++qidx) {
        search(index, queries[qidx], search_scheme, sm, internal_delegate);
    }
}

}
