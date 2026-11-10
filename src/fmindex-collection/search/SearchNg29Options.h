// SPDX-FileCopyrightText: 2026 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "CachedSearchScheme.h"
#include "Restore.h"
#include "SelectCursor.h"

#include <array>
#include <cstddef>
#include <map>
/**
 * recreating a new search, with an alignment matrix
 */
namespace fmc::search_ng29_options {

/**
 * Options:
 */
struct Options {
    bool indelCond{};
    bool matchCond{};
    bool reportCond{};
    bool duplicateDetection{};
    bool duplicateDetection2{};
};

class RangeSet {
public:
    std::map<size_t, size_t> ranges;
public:
    void clear() {
        ranges.clear();
    }
    void insert(size_t _pos, size_t _len) {
        auto iter = ranges.upper_bound(_pos+_len);
        while (iter != ranges.begin()) {
            auto prev = std::prev(iter);
            auto const& [p, l] = *prev;
            if (p+l < _pos) break;
            auto newPos = std::min(p, _pos);
            auto newLen = std::max(p+l, _pos + _len) - newPos;
            _pos = newPos;
            _len = newLen;
            ranges.erase(prev);
            iter = ranges.upper_bound(_pos + _len);
        }

        ranges[_pos] = _len;
    }

    bool contains(size_t _pos, size_t _len) const {
        auto iter = ranges.upper_bound(_pos+_len);
        if (iter == ranges.begin()) return false;
        auto const& [pos, len] = *std::prev(iter);
        return _pos >= pos && _pos + _len <= pos+len;
    }

    bool contains(size_t v) const {
        auto iter = ranges.upper_bound(v);
        if (iter == ranges.begin()) return false;
        auto const& [pos, len] = *std::prev(iter);
        return v < pos+len;
    }
};

struct BoundCheck {
    RangeSet knownHitsLeft;
    RangeSet knownHitsRight;
    std::unordered_map<size_t, size_t> knownHitsLeft2;
    std::unordered_map<size_t, size_t> knownHitsRight2;

    void clear() {
        knownHitsLeft.clear();
        knownHitsRight.clear();
    }

/*    template <typename cursor_t>
    void addHitImpl(cursor_t const& cur) {
        for (size_t i{0}; i < cur.len; ++i) {
            auto& lv = knownHitsLeft[cur.lb+i];
            auto& rv = knownHitsRight[cur.lbRev-i];
            lv = std::max(lv, cur.len-i);
            rv = std::max(rv, cur.len-i);
        }
    }*/

    template <typename cursor_t>
    void addHitLeftImpl(cursor_t const& cur, size_t errorsLeft) {
        knownHitsLeft.insert(cur.lb, cur.len);
/*        for (size_t i{0}; i < cur.len; ++i) {
            auto& lv = knownHitsLeft2[cur.lb+i];
            lv = std::max(lv, cur.len-i);
        }*/
        if (errorsLeft == 0) return;
        for (auto c : cur.extendLeft()) {
            addHitLeftImpl(c, errorsLeft-1);
        }
    }
    template <typename cursor_t>
    void addHitRightImpl(cursor_t const& cur, size_t errorsLeft) {
        knownHitsRight.insert(cur.lbRev, cur.len);
/*        for (size_t i{0}; i < cur.len; ++i) {
            auto& rv = knownHitsRight2[cur.lbRev+i];
            rv = std::max(rv, cur.len-i);
        }*/
        if (errorsLeft == 0) return;
        for (auto c : cur.extendRight()) {
            addHitRightImpl(c, errorsLeft-1);
        }
    }



    template <typename cursor_t>
    void addHit(cursor_t const& cur, size_t errorsLeft) {
        addHitLeftImpl(cur, errorsLeft);
        addHitRightImpl(cur, errorsLeft);
    }

    template <typename cursor_t>
    bool hasHitLeft(cursor_t const& cur) const {
        bool b1 = knownHitsLeft.contains(cur.lb, cur.len);
/*        bool b2 = [&]() {
            auto iter = knownHitsLeft2.find(cur.lb);
            if (iter == knownHitsLeft2.end()) return false;
            return cur.len <= iter->second;
        }();
        if (b1 != b2) std::print("hasHitLeft differs: {} {}\n", b1, b2);*/
        return b1;
    }

    template <typename cursor_t>
    bool hasHitRight(cursor_t const& cur) const {
        bool b1 =  knownHitsRight.contains(cur.lbRev, cur.len);
/*        bool b2 = [&]() {
            auto iter = knownHitsRight2.find(cur.lbRev);
            if (iter == knownHitsRight2.end()) return false;
            return cur.len <= iter->second;
        }();

        if (b1 != b2) {
            std::print("hasHitRight differs: {} {}\n", b1, b2);
            std::print("cur: {} {} {}\n", cur.lb, cur.lbRev, cur.len);
            auto iter = knownHitsRight.ranges.upper_bound(cur.lbRev +cur.len);
            if (iter == knownHitsRight.ranges.begin()) {
                std::print("no fitting pos\n");
            } else {
                {
                    auto const& [pos, len] = *iter;
                    fmt::print("info: {} {}\n", pos, len);
                }
                {
                    auto const& [pos, len] = *std::prev(iter);
                    fmt::print("info: {} {}\n", pos, len);
                }
                if (std::prev(iter) != knownHitsRight.ranges.begin())
                {
                    auto const& [pos, len] = *std::prev(std::prev(iter));
                    fmt::print("info: {} {}\n", pos, len);
                }

            }
        }*/
        return b1;
    }
};

template <typename query_t, typename search_t>
struct CursorFrontStorage {
    using T = uint8_t;
    std::vector<T> errors;
    query_t const& query;
    search_t const& search;
    T invalid = std::numeric_limits<T>::max();

    std::vector<bool> touchedLeft;
    std::vector<bool> touchedRight;

    std::vector<std::tuple<uint64_t, uint64_t>> extensionCount;

    CursorFrontStorage(query_t const& _query, search_t const& _search)
        : query{_query}
        , search{_search}
    {
        touchedLeft.resize(search.pi.size()+1);
        touchedRight.resize(search.pi.size()+1);
        for (size_t i{1}; i < search.pi.size()+1; ++i) {
            touchedLeft[i]  = touchedLeft[i-1]  || search.pi[i-1] == 0;
            touchedRight[i] = touchedRight[i-1] || search.pi[i-1] == search.pi.size()-1;
        }

        {
            uint64_t leftCt{};
            uint64_t rightCt{};
            extensionCount.emplace_back(leftCt, rightCt);

            for (size_t i{0}; i < query.size(); ++i) {
                if (isRight(i)) rightCt += 1;
                else leftCt += 1;
                extensionCount.emplace_back(leftCt, rightCt);
            }
        }
    }

    bool isRight(size_t i) const {
        bool r = (i == 0) || (search.pi[i-1] < search.pi[i]);
        return r;
    }


    void minSize(size_t i) {
        if (i < errors.size()) return;
        errors.resize(i, invalid);
    }
    auto operator[](size_t i) -> T& {
        return errors[i];
    }
    auto operator[](size_t i) const -> T {
        return errors[i];
    }
};

template <typename query_t, typename search_t>
struct CursorFront {
    using Storage = CursorFrontStorage<query_t, search_t>;
    using T = Storage::T;

    Storage& storage;
    int64_t offset{};
    size_t startI, endI;

/*    rhs.endI - rhs.offset = rhs.startI - offset

    1-0 == 0 - offset
    offset = rhs.startI - rhs.endI + rhs.offset
    offset  =  0        - 1        + 0*/

    CursorFront(CursorFront const& rhs)
        : storage{rhs.storage}
        , offset{rhs.offset - int64_t(rhs.range())}
        , startI{rhs.startI}
        , endI{rhs.endI}
    {
        storage.errors.resize(rhs.endI - rhs.offset);
        storage.minSize(endI - offset);
    }
    CursorFront(Storage& _storage)
        : storage{_storage}
    {
        storage.errors.reserve(100);
        storage.minSize(1);
        storage[0] = 0;
        startI = 0;
        endI = 1;
//        rolldown();
//        adjustStartAndEnd();
    }
    auto value(size_t i) -> auto& {
        return storage[i - offset];
    }
    auto value(size_t i) const -> auto const& {
        return storage[i - offset];
    }
    auto operator=(CursorFront const& rhs) -> CursorFront& {
        offset = rhs.offset;
        startI = rhs.startI;
        endI   = rhs.endI;
        return *this;
    }
    void adjustStartAndEnd() {
//        auto oldStartI = startI;
        // move end forwards
/*        while (endI-offset < storage.errors.size() && value(endI) != storage.invalid) {
            ++endI;
        }*/
        // move start forwards
        while (startI < endI && value(startI) == storage.invalid) {
            ++startI;
        }

    /*    if (oldStartI < startI) {
            auto delta = startI - oldStartI;
            for (size_t i{startI}; i < endI; ++i) {
                storage.errors[i-offset-delta] = storage.errors[i-offset];
            }
            offset = offset + delta;
        }*/
        storage.errors.resize(endI - offset);
    }

    void updateField(size_t i, T v) {
        if (!isInRange(i, v)) return;
        storage.minSize(i+1-offset);
        endI = std::max(endI, i+1);
        value(i) = std::min<T>(value(i), v);
    }

    bool isInRange(size_t i, T v) const {
        if (i == storage.query.size()) {
            if (v < storage.search.l.back()) return false;
        } else if (v < storage.search.l[i]) {
            return false;
        }

        if (i == 0) {
            if (v > 0) return false;
        } else if (v > storage.search.u[i-1]) {
            return false;
        }

        return true;
    }

    auto range() const -> size_t {
        return endI - startI;
    }

/*    void rolldown() {
        for (size_t i{startI}; i < endI; ++i) {
            if (value(i) == storage.invalid) continue;
            if (i == storage.query.size()) continue;
            updateField(i+1, value(i)+1);
        }
    }*/

    auto update(uint64_t symb, bool right) const -> CursorFront {
        auto cf = CursorFront{*this};

        // deletion
        for (size_t i{startI}; i < endI; ++i) {
            cf.value(i) = value(i)+1;
            if (!cf.isInRange(i, cf.value(i))) {
                cf.value(i) = cf.storage.invalid;
            }
        }
        // substitutions/matches
        for (size_t i{startI}; i < endI; ++i) {
            bool match = (symb == storage.query[i]);
            cf.updateField(i+1, value(i) + (match?0:1));
            if (isRight(i) != right) cf.value(i) = storage.invalid;
        }

        // insertions
        for (size_t i{cf.startI}; i < std::min(cf.endI, storage.query.size()); ++i) {
            if (cf.value(i) == storage.invalid) continue;
            cf.updateField(i+1, cf.value(i)+1);
        }

        cf.adjustStartAndEnd();
        return cf;
    }

    auto lowestValue() const -> T {
        T v = storage.invalid;
        for (size_t i{startI}; i < endI; ++i) {
            v = std::min(value(i), v);
        }
        return v;
    }
    bool isRight(size_t i) const {
        return storage.isRight(i);
    }
    auto requiresExtendLeftAndRight() const -> std::tuple<bool, bool> {
        auto [lct1, rct1] = storage.extensionCount[startI];
        auto [lct2, rct2] = storage.extensionCount[endI];
        return {lct1 != lct2, rct1 != rct2};
    }
};


template <bool Edit, typename index_t, typename query_t, typename search_t, typename delegate_t>
struct Search {
    constexpr static size_t Sigma = index_t::Sigma;
    constexpr static size_t FirstSymb = []() -> size_t {
        if constexpr (requires() { { index_t::FirstSymb }; }) {
            return index_t::FirstSymb;
        }
        return 1;
    }();

    using cursor_t = select_cursor_t<index_t>;
    using cursor_front_t = CursorFront<query_t, search_t>;

    index_t const& index;
    query_t const& query;
    search_t const& search;
    delegate_t const& delegate;

    Options options;

    BoundCheck& boundCheck;

    Search(index_t const& _index, query_t const& _query, search_t const& _search, delegate_t const& _delegate, Options _options, BoundCheck& _boundCheck)
        : index     {_index}
        , query     {_query}
        , search    {_search}
        , delegate  {_delegate}
        , options   {_options}
        , boundCheck {_boundCheck}
    {
    }

    bool run() {
        auto cur = cursor_t{index};
        auto storage = typename cursor_front_t::Storage{query, search};
        auto cf = cursor_front_t{storage};
        return search_next_symb(cur, cf);
    }

    bool search_next_symb(cursor_t const& cur, cursor_front_t& cf) const noexcept {
        if (cur.count() == 0) return false;
        auto r_cf_endI = Restore{cf.endI};

        if (options.duplicateDetection) {
/*            auto hitLeft = boundCheck.hasHitLeft(cur);
            auto hitRight = boundCheck.hasHitRight(cur);
            while (cf.endI > cf.startI) {
                auto touchLeft  = cf.storage.touchedLeft[cf.endI];
                auto touchRight = cf.storage.touchedRight[cf.endI];
                if (touchLeft && hitLeft) {
                    cf.endI -= 1;
                } else if (touchRight && hitRight) {
                    cf.endI -= 1;
                } else {
                    break;
                }
            }
            cf.storage.errors.resize(cf.endI-cf.offset);*/

            bool touchLeft  = cf.storage.touchedLeft[cf.startI];
            bool touchRight = cf.storage.touchedRight[cf.startI];

            if (touchLeft && boundCheck.hasHitLeft(cur)) {
                return false;
            }
            if (touchRight && boundCheck.hasHitRight(cur)) {
                return false;
            }
        }
        if (cf.range() == 0) return false;


        if (cf.endI == query.size()+1) {
            if (options.duplicateDetection) {
                if (boundCheck.hasHitLeft(cur) || boundCheck.hasHitRight(cur)) {
                    return false;
                }
                boundCheck.addHit(cur, 0);
//                boundCheck.addHit(cur, search.u.back() - cf.value(cf.endI-1));
            }

            auto e = cf.lowestValue();
            return delegate(cur, e);
        }

        // Optimization for cases the Cursor is only marking a single row
        if (cur.count() == 1) {
            auto [left, right] = cf.requiresExtendLeftAndRight();
            if (left) {
                auto [symb, newCur] = cur.extendLeftBySymbol();
                auto newCf = cf.update(symb, false);
                if (newCf.range() > 0) {
                    auto f = search_next_symb(newCur, newCf);
                    if (f) return true;
                }
            }
            if (right) {
                auto [symb, newCur] = cur.extendRightBySymbol();
                auto newCf = cf.update(symb, true);
                if (newCf.range() > 0) {
                    auto f = search_next_symb(newCur, newCf);
                    if (f) return true;
                }
            }
            return false;
        }
        // Optimization for the case, the alignment matrix only has a single valid entry
/*        if (cf.range() == 1) {
            if (cf.value(cf.startI) == search.u[cf.startI]) {
                auto symb = query[cf.startI];
                bool right = cf.isRight(cf.startI);
                auto newCur = right?cur.extendRight(symb):cur.extendLeft(symb);
                if (newCur.count() == 0) return false;
                auto newCf = cf.update(symb, right);
                if (newCf.range() == 0) return false;
                return search_next_symb(newCur, newCf);
            }
        }*/

        auto [left, right] = cf.requiresExtendLeftAndRight();

/*        auto next = [&](bool right) {
            for (size_t symb{FirstSymb}; symb < Sigma; ++symb) {
                auto newCur = [&]() {
                    if (right) return cur.extendRight(symb);
                    return cur.extendLeft(symb);
                }();
                if (newCur.count() == 0) continue;
                auto newCf = cf.update(symb, right);
                if (newCf.range() == 0) continue;
                auto f = search_next_symb(newCur, newCf);
                if (f) return true;
            }
            return false;
        };

        if (left) {
            auto f = next(false);
            if (f) return true;
        }
        if (right) {
            auto f = next(true);
            if (f) return true;
        }*/

        auto next = [&](auto const& cursors, bool right) {
            for (size_t symb{FirstSymb}; symb < Sigma; ++symb) {
                if (cursors[symb].count() == 0) continue;
                auto newCf = cf.update(symb, right);
                if (newCf.range() == 0) continue;
                auto f = search_next_symb(cursors[symb], newCf);
                if (f) return true;
            }
            return false;
        };

        if (left) {
            auto cursors = cur.extendLeft();
            auto f = next(cursors, false);
            if (f) return true;
        }
        if (right) {
            auto cursors = cur.extendRight();
            auto f = next(cursors, true);
            if (f) return true;
        }
        return false;
    }
};


template <bool Edit, typename index_t, Sequence query_t, typename delegate_t>
void search_impl(index_t const& index, query_t const& query, search_scheme::Scheme const& search_scheme, std::vector<size_t> const& partition, Options options, delegate_t&& delegate) {
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

    using TSearch = Search<Edit, index_t, query_t, decltype(*search_scheme.begin()), decltype(internal_delegate)>;

    // reorder query
    auto ordered_query = query;

    BoundCheck boundCheck;
    for (auto const& search : search_scheme) {
        for (size_t i{0}; i < query.size(); ++i) {
            ordered_query[i] = query[search.pi[i]];
        }
        bool f = TSearch{index, ordered_query, search, internal_delegate, options, boundCheck}.run();
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
void search_n_impl(index_t const& index, queries_t&& queries, selectSearchScheme_t&& selectSearchScheme, std::vector<size_t> const& partition, Options options, delegate_t&& delegate, size_t n) {
    if (queries.empty()) return;
    if (n == 0) return;
    for (size_t qidx{}; qidx < queries.size(); ++qidx) {
        size_t ct{};
        auto const& search_scheme = selectSearchScheme(queries[qidx].size());
        search_impl<Edit>(index, queries[qidx], search_scheme, partition, options, [&] (auto cur, size_t e) {
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
void search(index_t const& index, queries_t&& queries, search_scheme::Scheme const& search_scheme, std::vector<size_t> const& partition, Options options, delegate_t&& delegate, size_t n = std::numeric_limits<size_t>::max()) {
    // function that selects a search scheme
    auto selectSearchScheme = [&]([[maybe_unused]] size_t length) -> auto& {
        return search_scheme;
    };
    search_n_impl<Edit>(index, queries, selectSearchScheme, partition, options, delegate, n);
}

}
