#pragma once

#include "backtracking.h"
#include <cassert>

namespace oss::generator {

inline auto kianfar(int minK, int K) -> Scheme {
    assert(minK <= K);
    if (K == 0) {
        return {
            {{1}, {0}, {0}}
        };
    } else if (minK == 0 and K == 1) {
        return {
            {{1, 2}, {0, 0}, {0, 1}},
            {{2, 1}, {0, 0}, {0, 1}}
        };
    } else if (minK == 0 and K == 2) {
        return {
            {{1, 2, 3}, {0, 0, 2}, {0, 1, 2}},
            {{3, 2, 1}, {0, 0, 0}, {0, 2, 2}},
            {{2, 3, 1}, {0, 1, 1}, {0, 1, 2}}
        };
    } else if (minK == 0 and K == 3) {
        return {
            {{1, 2, 3, 4}, {0, 0, 0, 3}, {0, 2, 3, 3}},
            {{2, 3, 4, 1}, {0, 0, 0, 0}, {1, 2, 3, 3}},
            {{3, 4, 2, 1}, {0, 0, 2, 2}, {0, 0, 3, 3}},
        };
    }
    return {};
}

}
