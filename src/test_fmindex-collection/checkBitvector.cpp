// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include <catch2/catch_all.hpp>
#include <fmindex-collection/utils.h>

#include "allBitVectors.h"

TEMPLATE_TEST_CASE("check if rank on the bit vectors is working", "[BitVector]", ALLBITVECTORS) {
    using Vector = TestType;
    auto text = std::vector<uint8_t>{0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1};

    auto vec = Vector{std::span{text}};
    REQUIRE(vec.size() == text.size());

    SECTION("check that symbol() call works") {
        for (size_t i{0}; i < text.size(); ++i) {
            INFO(i);
            CHECK(vec.symbol(i) == text.at(i));
        }
    }

    SECTION("test complete vector on rank()") {
        CHECK(vec.rank( 0) == 0);
        CHECK(vec.rank( 1) == 0);
        CHECK(vec.rank( 2) == 1);
        CHECK(vec.rank( 3) == 2);
        CHECK(vec.rank( 4) == 2);
        CHECK(vec.rank( 5) == 2);
        CHECK(vec.rank( 6) == 3);
        CHECK(vec.rank( 7) == 3);
        CHECK(vec.rank( 8) == 4);
        CHECK(vec.rank( 9) == 5);
        CHECK(vec.rank(10) == 6);
        CHECK(vec.rank(11) == 6);
        CHECK(vec.rank(12) == 6);
        CHECK(vec.rank(13) == 6);
        CHECK(vec.rank(14) == 7);
    }
}
