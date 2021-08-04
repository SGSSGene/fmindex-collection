#include "suffixFilter.h"

#include <cassert>

namespace oss::generator {

auto suffixFilter(int N, int minK, int K) -> Scheme {
    assert(N > 0);
    assert(minK >= 0);
    assert(minK <= K);
    assert(K >= 0);
    assert(K < N);

    auto res = Scheme{};
    for (int n{0}; n < N; ++n) {
        auto s = SearchTree{};

        // generate suffix filter seeds
        for (int j{n}; j < N; ++j) {
            s.pi.push_back(j+1);
            s.l.push_back(0);
            s.u.push_back(std::min(j-n, K));

        }
        // fill rest of pattern
        for (int j{n-1}; j >= 0; --j) {
            s.pi.push_back(j+1);
            s.l.push_back(std::min(K, 1));
            s.u.push_back(K);
        }
        res.push_back(s);
    }

    // set minimum value
    for (auto& s : res) {
        s.l.back() = std::max(s.l.back(), minK);
    }

    return res;
}
}
