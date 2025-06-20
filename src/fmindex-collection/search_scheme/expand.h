// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "Scheme.h"
#include "isValid.h"
#include "nodeCount.h"
#include "weightedNodeCount.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <optional>
#include <type_traits>

namespace fmc::search_scheme {
namespace expand_detail {

inline auto forwards(std::vector<size_t> const& pi) -> std::vector<size_t> {
    auto rets = std::vector<size_t>{};
    rets.push_back(pi.size() == 1 or pi[1] > pi[0]?1:0);
    for (size_t i{1}; i < pi.size(); ++i) {
        rets.push_back(pi[i] > pi[i-1]?1:0);
    }
    return rets;
}
}

/**
 * Creates a partitioning of a search scheme of length oldLen to length newLen
 *
 * detail: returns a vector of length 'oldLen' which has `newLen` elements spread
 * uniformly across.
 */
inline auto expandCount(size_t oldLen, size_t newLen) -> std::vector<size_t> {
    assert(oldLen > 0);
    assert(newLen > 0);

    auto const block_size = newLen / oldLen;
    auto const rest       = newLen % oldLen;
    auto counts = std::vector<size_t>(oldLen, block_size);
    for (size_t i{0}; i < rest; ++i) {
        counts[i] += 1;
    }
    return counts;
}

inline bool isExpandable(Search s, size_t newLen) {
    assert(isValid(s));
    auto counts = expandCount(s.pi.size(), newLen);
    // check that all parts are still satisfiable

    if (s.l[0] > counts[s.pi[0]]) {
        return false;
    }
    for (size_t i{1}; i<s.pi.size(); ++i) {
        if (s.l[i]-s.u[i] > counts[s.pi[i]]) {
            return false;
        }
    }
    return true;
}


// Expands a pi vector according to 'counts'
inline auto expandPI(std::vector<size_t> const& pi, std::vector<size_t> const& counts) -> std::vector<size_t> {
    auto starts = std::vector<size_t>(pi.size(), 0);
    for (size_t i{1}; i < pi.size(); ++i) {
        starts[i] = starts[i-1] + counts[i-1];
    }
    auto nums = std::vector<size_t>{};
    auto expandForwards = [&](size_t l, size_t u) {
        for (size_t j = l; j <= u; ++j) {
            nums.push_back(j);
        }
    };
    auto expandBackwards = [&](size_t l, size_t u) {
        for (size_t j{u+1}; j > l; --j) {
            nums.push_back(j-1);
        }
    };
    auto expand = [&](size_t i, bool forward) {
        auto l = starts[pi[i]];
        auto u = l + counts[pi[i]]-1;
        if (forward) expandForwards(l, u);
        else expandBackwards(l, u);
    };

    auto fs = expand_detail::forwards(pi);

    for (size_t i{0}; i < pi.size(); ++i) {
        expand(i, fs[i]);
    }

    return nums;
}

inline auto expandPI(std::vector<size_t> const& pi, size_t _newLen) -> std::vector<size_t> {
    auto counts = expandCount(pi.size(), _newLen);
    return expandPI(pi, counts);
}

inline auto expandLowerBound(std::vector<size_t> const& pi, std::vector<size_t> const& bound, std::vector<size_t> const& counts) -> std::vector<size_t> {
    auto expandedBound = std::vector<size_t>{};

    for (size_t i{0}; i < pi.size(); ++i) {
        auto count = counts[pi[i]];
        while(count > 1) {
            --count;
            expandedBound.push_back((i>0)?bound[i-1]:0);
        }
        if (count > 0) {
            expandedBound.push_back(bound[i]);
        } else if (count == 0 and not expandedBound.empty()) {
            expandedBound.back() = bound[i];
        }
    }
    return expandedBound;
}

inline auto expandLowerBound(std::vector<size_t> const& pi, std::vector<size_t> const& bound, size_t _newLen) -> std::vector<size_t> {
    auto counts = expandCount(pi.size(), _newLen);
    return expandLowerBound(pi, bound, counts);
}

inline auto expandUpperBound(std::vector<size_t> const& pi, std::vector<size_t> const& bound, std::vector<size_t> const& counts) -> std::vector<size_t> {
    auto expandedBound = std::vector<size_t>{};

    for (size_t i{0}; i < pi.size(); ++i) {
        auto count = counts[pi[i]];
        while (count --> 0) {
            expandedBound.push_back(bound[i]);
        }
    }
    return expandedBound;
}


inline auto expandUpperBound(std::vector<size_t> const& pi, std::vector<size_t> const& bound, size_t _newLen) -> std::vector<size_t> {
    auto counts = expandCount(pi.size(), _newLen);
    return expandUpperBound(pi, bound, counts);
}

inline auto expand(Search const& s, size_t newLen) -> std::optional<Search> {
    auto r = Search{};
    r.pi = expandPI(s.pi, newLen);
    r.l  = expandLowerBound(s.pi, s.l, newLen);
    r.u  = expandUpperBound(s.pi, s.u, newLen);
    if (not isValid(r)) {
        return std::nullopt;
    }
    return {r};
}
inline auto expand(Scheme const& ss, size_t newLen) -> Scheme {
    auto r = Scheme{};
    for (auto const& s : ss) {
        auto o = expand(s, newLen);
        if (o) {
            r.push_back(*o);
        }
    }
    return r;
}

/** special expand function, that is given how every part should be expanded
 */
inline auto expand(Search const& s, std::vector<size_t> const& counts) -> std::optional<Search> {
    auto r = Search{};
    r.pi = expandPI(s.pi, counts);
    r.l  = expandLowerBound(s.pi, s.l, counts);
    r.u  = expandUpperBound(s.pi, s.u, counts);
    if (not isValid(r)) {
        return std::nullopt;
    }
    return {r};
}

inline auto expand(Scheme const& ss, std::vector<size_t> const& parts) -> Scheme {
    auto r = Scheme{};
    for (auto const& s : ss) {
        auto o = expand(s, parts);
        if (o) {
            r.push_back(*o);
        }
    }
    return r;
}

template <bool Edit=false>
auto expandByNC(Scheme const& ss, size_t _newLen, size_t sigma) -> Scheme {
    if (ss.size() == 0) return {};
    auto additionalPos = _newLen - ss[0].pi.size();
    auto counts = std::vector<size_t>(ss[0].pi.size(), 1);

    for (size_t i{0}; i<additionalPos; ++i) {
        double bestVal = std::numeric_limits<double>::max();
        size_t bestPos = 0;
        for (size_t j{0}; j < ss[0].pi.size(); ++j) {
            counts[j] += 1;
            auto ess = expand(ss, counts);
            counts[j] -= 1;
            auto f = nodeCount<Edit>(ess, sigma);
            if (f < bestVal) {
                bestVal = f;
                bestPos = j;
            }
        }
        counts[bestPos] += 1;
    }

    return expand(ss, counts);
}

/** Computes a 'count' vector according to best weighted node count
 */
template <bool Edit=false>
auto optimizeByWNC(Scheme const& ss, size_t _newLen, size_t sigma, size_t N) -> std::vector<size_t> {
    if (ss.size() == 0) return {};

    auto counts = std::vector<size_t>(ss[0].pi.size(), 1);
    auto additionalPos = _newLen - ss[0].pi.size();
    for (size_t i{0}; i<additionalPos; ++i) {
        double bestVal = std::numeric_limits<double>::max();
        size_t bestPos = 0;
        for (size_t j{0}; j < ss[0].pi.size(); ++j) {
            counts[j] += 1;
            auto ess = expand(ss, counts);
            counts[j] -= 1;
            auto f = weightedNodeCount<Edit>(ess, sigma, N);
            if (f < bestVal) {
                bestVal = f;
                bestPos = j;
            }
        }
        counts[bestPos] += 1;
    }
    return counts;
}

/** Expands from lower by adding position steps by step
 */
template <bool Edit=false>
auto expandByWNC(Scheme ss, size_t _newLen, size_t sigma, size_t N) -> Scheme {
    return expand(ss, optimizeByWNC(ss, _newLen, sigma, N));
}


template <bool Edit=false>
auto optimizeByWNCTopDown(Scheme const& _ss, size_t _newLen, size_t sigma, size_t N, size_t steps) -> std::vector<size_t> {
    if (_ss.size() == 0) return {};

    // uniform distribute part size
    auto counts = expandCount(_ss[0].pi.size(), _newLen);

    double lastVal = [&]() {
        auto ess = expand(_ss, counts);
        auto f = weightedNodeCount<Edit>(ess, sigma, N);
        return f;
    }();


    // loop through each part and try to improve it, by increasing it size
    while (true) {
        double bestVal = lastVal;
        size_t bestI1;
        size_t bestI2;
        for (size_t i1{0}; i1 < _ss[0].pi.size(); ++i1) {
            for (size_t i2{0}; i2 < _ss[0].pi.size(); ++i2) {
                if (i1 == i2) continue;
                if (counts[i1] <= steps) continue;
                counts[i1] -= steps;
                counts[i2] += steps;
                auto ess = expand(_ss, counts);
                auto f = weightedNodeCount<Edit>(ess, sigma, N);
                if (f < bestVal) {
                    bestI1 = i1;
                    bestI2 = i2;
                    bestVal = f;
                }
                counts[i1] += steps;
                counts[i2] -= steps;
            }
        }
        if (bestVal == lastVal) break;
        lastVal = bestVal;
        counts[bestI1] -= steps;
        counts[bestI2] += steps;
    }
    return counts;
}


template <bool Edit=false>
auto expandByWNCTopDown(Scheme const& ss, size_t _newLen, size_t sigma, size_t N, size_t steps) -> Scheme {
    return expand(ss, optimizeByWNCTopDown(ss, _newLen, sigma, N, steps));
}


inline auto limitToHamming(Search s) -> Search {
    auto len = s.pi.size();
    // limit can only be increased by one
    for (size_t i{len-1}; i>0; --i) {
        if (s.l[i] == 0) break;
        s.l[i-1] = std::max(s.l[i-1], s.l[i]-1);
    }
    // upper limit can only be increased by one
    for (size_t i{1}; i < len; ++i) {
        s.u[i] = std::min(s.u[i], s.u[i-1]+1);
    }
    return s;
}
inline auto limitToHamming(Scheme ss) -> Scheme {
    for (auto& s : ss) {
        s = limitToHamming(s);
    }
    return ss;
}

}
