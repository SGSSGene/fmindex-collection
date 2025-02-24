// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include <catch2/catch_all.hpp>
#include <fmindex-collection/bitset_popcount.h>

#include <fstream>
#include <nanobench.h>

namespace fmc = fmindex_collection;

TEST_CASE("check if signed_rshift_and_count works", "[signed_rshift_and_count]") {
    SECTION("all ones") {
        auto b = std::bitset<64>{};
        b.flip();
        for (size_t i{0}; i < 65; ++i) {
            CHECK(fmc::signed_rshift_and_count(b, i) == 64-i);
        }

        for (size_t i{64}; i < 128; ++i) {
            CHECK(fmc::signed_rshift_and_count(b, i) == i-64);
        }
    }
}
TEST_CASE("count_first_or_last_n_bits", "[count_bits]") {
    SECTION("all ones") {
        auto b = std::bitset<64>{};
        b.flip();
        for (size_t i{0}; i < 65; ++i) {
            CHECK(fmc::skip_first_or_last_n_bits_and_count(b, i) == 64-i);
        }

        for (size_t i{64}; i <= 128; ++i) {
            CHECK(fmc::skip_first_or_last_n_bits_and_count(b, i) == i-64);
        }
    }
    SECTION("single one on the first bit") {
        auto b = std::bitset<64>{};
        b[0] = true;
        CHECK(fmc::skip_first_or_last_n_bits_and_count(b, 0) == 1);
        for (size_t i{1}; i <= 64; ++i) {
            CHECK(fmc::skip_first_or_last_n_bits_and_count(b, i) == 0);
        }
        for (size_t i{65}; i <= 128; ++i) {
            CHECK(fmc::skip_first_or_last_n_bits_and_count(b, i) == 1);
        }
    }

    SECTION("single one on the last bit") {
        auto b = std::bitset<64>{};
        b[63] = true;
        for (size_t i{0}; i < 64; ++i) {
            CHECK(fmc::skip_first_or_last_n_bits_and_count(b, i) == 1);
        }
        for (size_t i{64}; i < 128; ++i) {
            CHECK(fmc::skip_first_or_last_n_bits_and_count(b, i) == 0);
        }
        CHECK(fmc::skip_first_or_last_n_bits_and_count(b, 128) == 1);
    }
}
