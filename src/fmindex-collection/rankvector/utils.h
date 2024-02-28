// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>

namespace fmindex_collection::rankvector {

/** Counts how many bits are needed to represent the number y
 * It performs the computation `ceil(log_2(x))` if x > 0
 *                             `1` if x == 0
 */
constexpr inline uint64_t required_bits(uint64_t x) {
    if (x == 0) return 1;
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

static_assert(required_bits(16) == 5);
static_assert(required_bits(15) == 4);
static_assert(required_bits(8) == 4);
static_assert(required_bits(7) == 3);
static_assert(required_bits(4) == 3);
static_assert(required_bits(3) == 2);
static_assert(required_bits(2) == 2);
static_assert(required_bits(1) == 1);
static_assert(required_bits(0) == 1);

}
