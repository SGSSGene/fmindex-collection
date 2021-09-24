#pragma once

#include "../Scheme.h"
#include <cassert>

namespace search_schemes::generator {

inline auto kucherov(int N, int K) -> Scheme {
    assert((N == K+1) or (N == K+2));

    if (K == 0) {
        return {
            {{0}, {0}, {0}}
        };
    } else if (N == 2 and K == 1) {
        return {
            {{0, 1}, {0, 0}, {0, 1}},
            {{1, 0}, {0, 0}, {0, 1}}
        };
    } else if (N == 3 and K == 2) {
        return {
            {{0, 1, 2}, {0, 0, 0}, {0, 2, 2}},
            {{2, 1, 0}, {0, 0, 0}, {0, 1, 2}},
            {{1, 0, 2}, {0, 0, 1}, {0, 1, 2}}
        };
    } else if (N == 4 and K == 3) {
        return {
            {{0, 1, 2, 3}, {0, 0, 0, 0}, {0, 1, 3, 3}},
            {{1, 0, 2, 3}, {0, 0, 1, 1}, {0, 1, 3, 3}},
            {{2, 3, 1, 0}, {0, 0, 0, 0}, {0, 1, 3, 3}},
            {{3, 2, 1, 0}, {0, 0, 1, 1}, {0, 1, 3, 3}},
        };
    } else if (N == 5 and K == 4) {
        return {
            {{0, 1, 2, 3, 4}, {0, 0, 0, 0, 0}, {0, 2, 2, 4, 4}},
            {{4, 3, 2, 1, 0}, {0, 0, 0, 0, 0}, {0, 1, 3, 4, 4}},
            {{1, 0, 2, 3, 4}, {0, 0, 1, 3, 3}, {0, 1, 3, 4, 4}},
            {{0, 1, 2, 3, 4}, {0, 0, 1, 3, 3}, {0, 1, 3, 4, 4}},
            {{3, 2, 4, 1, 0}, {0, 0, 0, 1, 1}, {0, 1, 2, 4, 4}},
            {{2, 1, 0, 3, 4}, {0, 0, 0, 1, 3}, {0, 1, 2, 4, 4}},
            {{1, 0, 2, 3, 4}, {0, 0, 1, 2, 4}, {0, 1, 2, 4, 4}},
            {{0, 1, 2, 3, 4}, {0, 0, 0, 3, 4}, {0, 0, 4, 4, 4}},

        };
    } else if (N == 2 and K == 0) {
        return {
            {{0, 1}, {0, 0}, {0, 0}}
        };
    } else if (N == 3 and K == 1) {
        return {
            {{0, 1, 2}, {0, 0, 0}, {0, 1, 1}},
            {{1, 2, 0}, {0, 0, 0}, {0, 0, 1}}
        };
    } else if (N == 4 and K == 2) {
        return {
            {{0, 1, 2, 3}, {0, 0, 0, 0}, {0, 1, 1, 2}},
            {{3, 2, 1, 0}, {0, 0, 0, 0}, {0, 1, 2, 2}},
            {{1, 2, 3, 0}, {0, 0, 0, 1}, {0, 0, 1, 2}},
            {{0, 1, 2, 3}, {0, 0, 0, 2}, {0, 0, 2, 2}}
        };
    } else if (N == 5 and K == 3) {
        return {
            {{0, 1, 2, 3, 4}, {0, 0, 0, 0, 0}, {0, 1, 2, 3, 3}},
            {{1, 2, 3, 4, 0}, {0, 0, 0, 0, 0}, {0, 1, 2, 2, 3}},
            {{2, 3, 4, 1, 0}, {0, 0, 0, 0, 1}, {0, 1, 1, 3, 3}},
            {{3, 4, 2, 1, 0}, {0, 0, 0, 1, 2}, {0, 0, 3, 3, 3}},
        };
    } else if (N == 6 and K == 4) {
        return {
            {{0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 0, 0}, {0, 1, 2, 3, 4, 4}},
            {{1, 2, 3, 4, 5, 0}, {0, 0, 0, 0, 0, 0}, {0, 1, 2, 3, 4, 4}},
            {{5, 4, 3, 2, 1, 0}, {0, 0, 0, 0, 0, 1}, {0, 1, 2, 2, 4, 4}},
            {{3, 4, 5, 2, 1, 0}, {0, 0, 0, 0, 1, 2}, {0, 1, 1, 3, 4, 4}},
            {{2, 3, 4, 5, 1, 0}, {0, 0, 0, 0, 2, 3}, {0, 1, 1, 2, 4, 4}},
            {{4, 5, 3, 2, 1, 0}, {0, 0, 0, 1, 3, 3}, {0, 0, 3, 3, 4, 4}},
            {{0, 1, 2, 3, 4, 5}, {0, 0, 0, 3, 3, 3}, {0, 0, 3, 3, 4, 4}},
            {{0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 4, 4}, {0, 0, 2, 4, 4, 4}},
            {{2, 3, 1, 0, 4, 5}, {0, 0, 0, 1, 2, 4}, {0, 0, 2, 2, 4, 4}},
            {{4, 5, 3, 2, 1, 0}, {0, 0, 0, 0, 4, 4}, {0, 0, 1, 4, 4, 4}},
        };
    }

    return {};
}

}