// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>

namespace fmindex_collection::occtable {

/** Counts how many bits are needed to represent the number y
 * It performs the computation `ceil(log_2(x))` if x > 0
 *                             `0` if x == 0
 */
constexpr inline uint64_t required_bits(uint64_t x) {
    uint64_t i{0};
    while (x != 0) {
        x = x >> 1;
        ++i;
    }
    return i;
}

// computes b to the power of n
constexpr inline uint64_t pow(uint64_t b, uint64_t n) {
    if (n == 0) return 1;
    return pow(b, (n-1)) * b;
}

}
