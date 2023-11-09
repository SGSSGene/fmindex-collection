// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "Scheme.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>

namespace search_schemes {

/**
 * \tparam Edit use edit distance, other wise Hamming distance
 * \param s search scheme
 * \param sigma size of the alphabet (without delimiter)
 * \param N     size of the reference text, ~3'000'000'000 for hg
 */
template <bool Edit>
long double weightedNodeCount(Search s, size_t sigma, size_t N) {
    auto n_max = s.pi.size();
    auto e     = *std::max_element(begin(s.u), end(s.u));

    auto lastArray = std::vector<long double>(e+1, 0);
    lastArray[0] = 1;

    long double acc = 0;

    auto newArray = std::vector<long double>(e+1, 0);
    for (size_t n {1}; n <= n_max; ++n) {
        auto f = N / std::pow(sigma, n);
        if (f > 1) f = 1.;

        for (size_t i{0}; i < e+1; ++i) {
            if (s.l[n-1] <= i and i <= s.u[n-1]) {
                newArray[i] = lastArray[i];
                if (i > 0) {
                    if constexpr (Edit) {
                        newArray[i] += (sigma-1) * lastArray[i-1] + (sigma) * lastArray[i-1] + lastArray[i-1];
                    } else {
                        newArray[i] += (sigma-1) * lastArray[i-1];
                    }
                }
                newArray[i] *= f;
                acc += newArray[i];

            } else {
                newArray[i] = 0;
            }
        }
        std::swap(newArray, lastArray);
    }
    return acc;
}

/**
 * \tparam Edit use edit distance, other wise Hamming distance
 * \param ss search schemes
 * \param sigma size of the alphabet (without delimiter)
 * \param N     size of the reference text, ~3'000'000'000 for hg
 */
template <bool Edit>
long double weightedNodeCount(Scheme const& ss, size_t sigma, size_t N) {
    return std::accumulate(begin(ss), end(ss), static_cast<long double>(0.), [&](long double v, auto const& s) {
        return v + weightedNodeCount<Edit>(s, sigma, N);
    });
}
}
