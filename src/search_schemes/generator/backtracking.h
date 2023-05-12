// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
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
