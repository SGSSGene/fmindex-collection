#pragma once

#include "../Scheme.h"
#include <cassert>

namespace search_schemes::generator {

inline auto kianfar(int K) -> Scheme {
    if (K == 0) {
        return {
            {{0}, {0}, {0}}
        };
    } else if (K == 1) {
        return {
            {{0, 1}, {0, 0}, {0, 1}},
            {{1, 0}, {0, 1}, {0, 1}}
        };
    } else if (K == 2) {
        return {
            {{0, 1, 2}, {0, 0, 2}, {0, 1, 2}},
            {{2, 1, 0}, {0, 0, 0}, {0, 2, 2}},
            {{1, 2, 0}, {0, 1, 1}, {0, 1, 2}}
        };
    } else if (K == 3) {
        return {
            {{0, 1, 2, 3}, {0, 0, 0, 3}, {0, 2, 3, 3}},
            {{1, 2, 3, 0}, {0, 0, 0, 0}, {1, 2, 3, 3}},
            {{2, 3, 1, 0}, {0, 0, 2, 2}, {0, 0, 3, 3}},
        };
    }
    return {};
}

}
