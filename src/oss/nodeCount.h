#pragma once

#include "Scheme.h"

#include <algorithm>
#include <numeric>

namespace oss {

inline long double nodeCount(SearchTree s, int sigma) {
    auto n_max = s.pi.size();
    auto e     = *std::max_element(begin(s.u), end(s.u));

    auto lastArray = std::vector<long double>(e+1, 0);
    lastArray[0] = 1;

    long double acc = 0;

    auto newArray = std::vector<long double>(e+1, 0);
    for (int n {1}; n <= n_max; ++n) {
        for (int i{0}; i < e+1; ++i) {
            if (s.l[n-1] <= i and i <= s.u[n-1]) {
                newArray[i] = lastArray[i];
                if (i > 0) {
                    newArray[i] += (sigma-1) * lastArray[i-1];
                }
                acc += newArray[i];
            } else {
                newArray[i] = 0;
            }
        }
        std::swap(newArray, lastArray);
    }
    return acc;
}

inline long double nodeCount(Scheme const& ss, int sigma) {
    return std::accumulate(begin(ss), end(ss), (long double){0.}, [&](long double v, auto const& s) {
        return v + nodeCount(s, sigma);
    });
}

inline long double nodeCountEdit(SearchTree s, int sigma) {
    auto n_max = s.pi.size();
    auto e     = *std::max_element(begin(s.u), end(s.u));

    auto lastArray = std::vector<long double>(e+1, 0);
    lastArray[0] = 1;

    long double acc = 0;

    auto newArray = std::vector<long double>(e+1, 0);
    for (int n {1}; n <= n_max; ++n) {
        for (int i{0}; i < e+1; ++i) {
            if (s.l[n-1] <= i and i <= s.u[n-1]) {
                newArray[i] = lastArray[i];
                if (i > 0) {
                    newArray[i] += (sigma-1) * lastArray[i-1] + (sigma) * lastArray[i-1] + lastArray[i-1];
                }
                acc += newArray[i];
            } else {
                newArray[i] = 0;
            }
        }
        std::swap(newArray, lastArray);
    }
    return acc;
}

inline long double nodeCountEdit(Scheme const& ss, int sigma) {
    return std::accumulate(begin(ss), end(ss), (long double){0.}, [&](long double v, auto const& s) {
        return v + nodeCountEdit(s, sigma);
    });
}

}
