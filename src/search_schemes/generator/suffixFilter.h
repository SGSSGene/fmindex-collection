// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../Scheme.h"

#include <cassert>

namespace search_schemes::generator {

inline auto suffixFilter(size_t N, size_t minK, size_t K) -> Scheme {
    assert(N > 0);
    assert(minK <= K);
    assert(K < N);

    auto res = Scheme{};
    for (size_t n{0}; n < N; ++n) {
        auto s = Search{};

        // generate suffix filter seeds
        for (size_t j{n}; j < N; ++j) {
            s.pi.push_back(j);
            s.l.push_back(0);
            s.u.push_back(std::min(j-n, K));

        }
        // fill rest of pattern
        for (size_t j{n}; j > 0; --j) {
            s.pi.push_back(j-1);
            s.l.push_back(std::min(K, size_t{1}));
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
