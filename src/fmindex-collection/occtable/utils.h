// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

namespace fmindex_collection {
namespace occtable {

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
}
