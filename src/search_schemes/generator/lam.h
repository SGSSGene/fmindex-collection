// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../Scheme.h"
#include <cassert>

namespace search_schemes::generator {

inline auto lam(int K) -> Scheme {
    assert(K==2);
    if (K == 2) {
        return {
            {{0, 1, 2}, {0, 0, 0}, {0, 2, 2}},
            {{2, 1, 0}, {0, 0, 0}, {0, 1, 2}},
            {{1, 2, 0}, {0, 0, 1}, {0, 1, 2}}
        };
    }

    return {};
}

}
