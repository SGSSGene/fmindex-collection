// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../Scheme.h"

#include <cassert>
#include <numeric>


namespace search_schemes::generator {

inline auto backtracking(size_t N, size_t minK, size_t K) -> Scheme {
    assert(N > 0);
    assert(K >= minK);
    auto s = Search{std::vector<size_t>(N, 0), std::vector<size_t>(N, 0), std::vector<size_t>(N, K)};
    s.l.back() = minK;
    std::iota(begin(s.pi), end(s.pi), 0);
    return Scheme{s};
}

}
