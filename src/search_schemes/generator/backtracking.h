#pragma once

#include "../Scheme.h"
#include <numeric>


namespace search_schemes::generator {

inline auto backtracking(int N, int minK, int K) -> Scheme {
    auto s = Search{std::vector<int>(N, 0), std::vector<int>(N, 0), std::vector<int>(N, K)};
    s.l.back() = minK;
    std::iota(begin(s.pi), end(s.pi), 1);
    return Scheme{s};
}

}
