// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "Scheme.h"

#include <algorithm>
#include <cassert>
#include <numeric>

namespace search_schemes {

/**
 * If its not edit distance, its Hamming distance
 */
template <bool Edit>
long double nodeCount(Search s, size_t sigma) {
    auto n_max = s.pi.size();
    auto e     = *std::max_element(begin(s.u), end(s.u));

    auto lastArray = std::vector<long double>(e+1, 0);
    lastArray[0] = 1;

    long double acc = 0;

    auto newArray = std::vector<long double>(e+1, 0);
    for (size_t n {1}; n <= n_max; ++n) {
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
 * If its not edit distance, its Hamming distance
 */
template <bool Edit>
long double nodeCount(Scheme const& ss, size_t sigma) {
    return std::accumulate(begin(ss), end(ss), static_cast<long double>(0.), [&](long double v, auto const& s) {
        return v + nodeCount<Edit>(s, sigma);
    });
}
}
