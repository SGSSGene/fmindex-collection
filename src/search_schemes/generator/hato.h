// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../Scheme.h"
#include <cassert>

namespace search_schemes::generator {

inline auto hato(int K) -> Scheme {
    if (K == 0) {
        return {
            {{0}, {0}, {0}} // not found by hato
        };
    } else if (K == 1) { // not found by hatto
        return {
            {{0, 1}, {0, 0}, {0, 1}},
            {{1, 0}, {0, 0}, {0, 1}}
        };
    } else if (K == 2) {
        return {
            {{0, 1, 2}, {0, 1, 1}, {0, 2, 2}},
            {{1, 0, 2}, {0, 0, 0}, {0, 1, 2}},
            {{2, 1, 0}, {0, 0, 2}, {0, 1, 2}}
        };
    } else if (K == 3) {
        return {
            {{0, 1, 2, 3}, {0, 0, 0, 0}, {0, 1, 3, 3}},
            {{1, 0, 2, 3}, {0, 1, 1, 1}, {0, 1, 3, 3}},
            {{2, 3, 1, 0}, {0, 0, 0, 2}, {0, 1, 3, 3}},
            {{3, 2, 1, 0}, {0, 1, 1, 3}, {0, 1, 3, 3}},
        };
    } else if (K == 4) {
        return {
            {{0, 1, 2, 3, 4}, {0, 0, 2, 2, 2}, {0, 2, 2, 4, 4}},
            {{1, 2, 0, 3, 4}, {0, 0, 0, 0, 0}, {0, 1, 2, 4, 4}},
            {{2, 1, 0, 3, 4}, {0, 1, 1, 1, 1}, {0, 1, 2, 4, 4}},
            {{3, 4, 2, 1, 0}, {0, 0, 0, 0, 3}, {0, 1, 4, 4, 4}},
            {{4, 3, 2, 1, 0}, {0, 1, 1, 1, 4}, {0, 1, 4, 4, 4}},
        };
    } else if (K == 5) {
        return {
            {{0, 1, 2, 3, 4, 5}, {0, 0, 0, 2, 2, 2}, {0, 1, 3, 5, 5, 5}},
            {{1, 0, 2, 3, 4, 5}, {0, 1, 1, 3, 3, 3}, {0, 1, 3, 5, 5, 5}},
            {{2, 3, 1, 0, 4, 5}, {0, 0, 0, 0, 0, 0}, {0, 1, 3, 3, 5, 5}},
            {{3, 2, 1, 0, 4, 5}, {0, 1, 1, 1, 1, 1}, {0, 1, 3, 3, 5, 5}},
            {{4, 5, 3, 2, 1, 0}, {0, 0, 0, 0, 0, 4}, {0, 1, 3, 5, 5, 5}},
            {{5, 4, 3, 2, 1, 0}, {0, 1, 1, 1, 1, 5}, {0, 1, 3, 5, 5, 5}},
        };
    } else if (K == 6) {
        return {
            {{0, 1, 2, 3, 4, 5, 6}, {0, 0, 2, 2, 2, 2, 6}, {0, 2, 2, 6, 6, 6, 6}},
            {{1, 2, 0, 3, 4, 5, 6}, {0, 1, 1, 1, 1, 1, 5}, {0, 1, 2, 6, 6, 6, 6}},
            {{2, 1, 0, 3, 4, 5, 6}, {0, 0, 0, 0, 0, 0, 4}, {0, 1, 2, 6, 6, 6, 6}},
            {{3, 4, 5, 6, 2, 1, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 1, 3, 3, 6, 6, 6}},
            {{4, 3, 5, 6, 2, 1, 0}, {0, 1, 1, 1, 1, 1, 1}, {0, 1, 3, 3, 6, 6, 6}},
            {{5, 6, 4, 3, 2, 1, 0}, {0, 0, 0, 2, 2, 2, 2}, {0, 1, 3, 3, 6, 6, 6}},
            {{6, 5, 4, 3, 2, 1, 0}, {0, 1, 1, 3, 3, 3, 3}, {0, 1, 3, 3, 6, 6, 6}},
        };
    } else if (K == 7) {
        return {
            {{0, 1, 2, 3, 4, 5, 6, 7}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 3, 3, 7, 7, 7, 7}},
            {{1, 0, 2, 3, 4, 5, 6, 7}, {0, 1, 1, 1, 1, 1, 1, 1}, {0, 1, 3, 3, 7, 7, 7, 7}},
            {{2, 3, 1, 0, 4, 5, 6, 7}, {0, 0, 0, 2, 2, 2, 2, 2}, {0, 1, 3, 3, 7, 7, 7, 7}},
            {{3, 2, 1, 0, 4, 5, 6, 7}, {0, 1, 1, 3, 3, 3, 3, 3}, {0, 1, 3, 3, 7, 7, 7, 7}},
            {{4, 5, 6, 7, 3, 2, 1, 0}, {0, 0, 0, 0, 0, 0, 0, 4}, {0, 1, 3, 3, 7, 7, 7, 7}},
            {{5, 4, 6, 7, 3, 2, 1, 0}, {0, 1, 1, 1, 1, 1, 1, 5}, {0, 1, 3, 3, 7, 7, 7, 7}},
            {{6, 7, 5, 4, 3, 2, 1, 0}, {0, 0, 0, 2, 2, 2, 2, 6}, {0, 1, 3, 3, 7, 7, 7, 7}},
            {{7, 6, 5, 4, 3, 2, 1, 0}, {0, 1, 1, 3, 3, 3, 3, 7}, {0, 1, 3, 3, 7, 7, 7, 7}},
        };
    }

    return {};
}

}