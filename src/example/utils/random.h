// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <cstdint>

static uint32_t x=123456789, y=362436069, z=521288629;

inline void xorshf96_reset() {
    x=123456789;
    y=362436069;
    z=521288629;
}
inline uint32_t xorshf96() {          //period 2^96-1
    uint32_t t;
    x ^= x << 16;
    x ^= x >> 5;
    x ^= x << 1;

   t = x;
   x = y;
   y = z;
   z = t ^ x ^ y;

  return z;
}

