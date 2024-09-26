// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include <catch2/catch_all.hpp>
#include <fmindex-collection/bitset_popcount.h>

#include <fstream>
#include <nanobench.h>

namespace fmc = fmindex_collection;

TEST_CASE("check if lshift_and_count works", "[leftshift_and_count]") {

    auto b = std::bitset<64>{};
    b.flip();
    for (size_t i{0}; i < 65; ++i) {
        CHECK(fmc::signed_rshift_and_count(b, i) == 64-i);
    }

    for (size_t i{65}; i < 128; ++i) {
        CHECK(fmc::signed_rshift_and_count(b, i) == i-64);
    }
}
