// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "backtracking.h"
#include <cassert>

namespace fmc::search_scheme::generator {

inline auto bestKnown(size_t N, size_t minK, size_t K) -> Scheme {
    assert(N > 0);
    assert(minK <= K);
    if (N == 1 and K == 0) {
        return {
            {{0}, {0}, {0}}
        };
    } else if (N == 2 and minK == 0 and K == 1) {
        return {
            {{0, 1}, {0, 0}, {0, 1}},
            {{1, 0}, {0, 1}, {0, 1}}
        };
    } else if (N == 2 and minK == 1 and K == 1) {
        return {
            {{0, 1}, {0, 1}, {0, 1}},
            {{1, 0}, {0, 1}, {0, 1}}
        };
    } else if (N == 3 and minK == 0 and K == 2) {
        return {
            {{0, 1, 2, 3}, {0, 0, 1, 1}, {0, 0, 2, 2}},
            {{2, 1, 0, 3}, {0, 0, 0, 0}, {0, 1, 1, 2}},
            {{3, 2, 1, 0}, {0, 0, 0, 2}, {0, 1, 2, 2}}
        };
    } else if (N == 3 and minK == 1 and K == 2) {
        return {
            {{0, 1, 2, 3}, {0, 0, 0, 1}, {0, 0, 2, 2}},
            {{2, 1, 0, 3}, {0, 0, 1, 1}, {0, 1, 1, 2}},
            {{3, 2, 1, 0}, {0, 0, 0, 2}, {0, 1, 2, 2}}
        };
    } else if (N == 3 and minK == 2 and K == 2) {
        return {
            {{3, 2, 1, 0}, {0, 0, 1, 2}, {0, 0, 2, 2}},
            {{1, 2, 3, 0}, {0, 0, 0, 2}, {0, 1, 1, 2}},
            {{0, 1, 2, 3}, {0, 0, 0, 2}, {0, 1, 2, 2}}
        };
    } else if (N == 4 and minK == 0 and K == 3) {
        return {
            {{4, 3, 2, 1, 0}, {0, 0, 0, 0, 0}, {0, 0, 3, 3, 3}},
            {{2, 3, 4, 1, 0}, {0, 0, 1, 1, 1}, {0, 1, 1, 2, 3}},
            {{1, 2, 3, 4, 0}, {0, 0, 0, 2, 2}, {0, 1, 2, 2, 3}},
            {{0, 1, 2, 3, 4}, {0, 0, 0, 0, 3}, {0, 2, 2, 3, 3}}
        };
    } else if (N == 4 and minK == 1 and K == 3) {
        return {
            {{4, 3, 2, 1, 0}, {0, 0, 0, 0, 1}, {0, 0, 3, 3, 3}},
            {{2, 3, 4, 1, 0}, {0, 0, 1, 1, 1}, {0, 1, 1, 2, 3}},
            {{1, 2, 3, 4, 0}, {0, 0, 0, 2, 2}, {0, 1, 2, 2, 3}},
            {{0, 1, 2, 3, 4}, {0, 0, 0, 0, 3}, {0, 2, 2, 3, 3}}
        };
    } else if (N == 4 and minK == 2 and K == 3) {
        return {
            {{4, 3, 2, 1, 0}, {0, 0, 0, 0, 2}, {0, 0, 3, 3, 3}},
            {{2, 3, 4, 1, 0}, {0, 0, 1, 1, 2}, {0, 1, 1, 2, 3}},
            {{1, 2, 3, 4, 0}, {0, 0, 0, 2, 2}, {0, 1, 2, 2, 3}},
            {{0, 1, 2, 3, 4}, {0, 0, 0, 0, 3}, {0, 2, 2, 3, 3}}
        };
    } else if (N == 4 and minK == 3 and K == 3) {
        return {
            {{4, 3, 2, 1, 0}, {0, 0, 0, 0, 3}, {0, 0, 3, 3, 3}},
            {{2, 3, 4, 1, 0}, {0, 0, 1, 1, 3}, {0, 1, 1, 2, 3}},
            {{1, 2, 3, 4, 0}, {0, 0, 0, 2, 3}, {0, 1, 2, 2, 3}},
            {{0, 1, 2, 3, 4}, {0, 0, 0, 0, 3}, {0, 2, 2, 3, 3}}
        };
    } else if (N == 5 and K == 4) {
        return {
            {{0, 1, 2, 3, 4}, {0, 0, 0, 0, std::max(minK, size_t{4})}, {0, 3, 3, 4, 4}},
            {{1, 2, 3, 4, 0}, {0, 0, 0, 0, std::max(minK, size_t{0})}, {2, 2, 3, 3, 4}},
            {{4, 3, 2, 1, 0}, {0, 0, 0, 3, std::max(minK, size_t{3})}, {0, 0, 4, 4, 4}}
        };
    }
    return backtracking(N, minK, K);
}

}
