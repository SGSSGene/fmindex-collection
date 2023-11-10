// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "Scheme.h"

#include <algorithm>
#include <cassert>
#include <cstdint>

namespace search_schemes {

namespace complete::detail {

using ErrorConfig = std::vector<uint8_t>;

inline auto covers(Search const& s, ErrorConfig const& errorConfig) -> bool {
    assert(s.pi.size() == errorConfig.size());
    size_t a{};
    for (size_t i{0}; i < s.pi.size(); ++i) {
        a += errorConfig.at(s.pi[i]);
        if (!(s.l[i] <= a && a <= s.u[i])) {
            return false;
        }
    }
    return true;
}

inline auto covers(Scheme const& ss, ErrorConfig const& errorConfig) -> bool {
    for (auto const& s : ss) {
        if (covers(s, errorConfig)) {
            return true;
        }
    }
    return false;
}

template <typename CB>
void generateErrorConfig(CB cb, ErrorConfig& errorConfig, size_t currentK, size_t startI, size_t maxK) {
    if (currentK >= maxK) return;
    for (size_t i{startI}; i < errorConfig.size(); ++i) {
        errorConfig[i] += 1;
        cb(errorConfig, currentK+1);
        generateErrorConfig(cb, errorConfig, currentK+1, i, maxK);
        errorConfig[i] -= 1;
    }
}


template <typename CB>
void generateErrorConfig(CB cb, size_t len, size_t minK, size_t maxK) {
    ErrorConfig errorConfig(len, 0);

    if (minK == 0) {
        cb(errorConfig);
    }
    generateErrorConfig([&](ErrorConfig const& config, size_t k) {
        if (k >= minK) {
            cb(config);
        }
    }, errorConfig, 0, 0, maxK);
}


/* checks if Scheme covers patterns up to a certain error
 *
 */
inline auto isComplete(Scheme const& ss, size_t minK, size_t maxK) -> bool {
    if (ss.empty()) return false;
    bool complete{true};
    auto len = ss.at(0).pi.size();
    generateErrorConfig([&](ErrorConfig const& config) {
        bool r = covers(ss, config);
        complete = complete && r;
    }, len, minK, maxK);
    return complete;
}

}

inline auto isComplete(Scheme const& ss, size_t minK, size_t maxK) -> bool {
    return complete::detail::isComplete(ss, minK, maxK);
}


}
