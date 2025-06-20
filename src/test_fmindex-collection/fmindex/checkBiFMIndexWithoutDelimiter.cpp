// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "../string/utils.h"

#include <algorithm>
#include <catch2/catch_all.hpp>
#include <fmindex-collection/fmindex/BiFMIndex.h>
#include <fmindex-collection/fmindex/BiFMIndexCursor.h>
#include <fstream>

TEST_CASE("checking bidirectional fm index without delimiters", "[bifmindex-nd]") {
    auto input = std::vector<std::vector<uint8_t>> {std::vector<uint8_t>{1, 2, 0}};
    SECTION("test index without delimiter") {
        auto index = fmindex_collection::BiFMIndex<255>{input, /*.samplingRate=*/ 1, /*.threadNbr=*/1, false};

        auto cursor = fmindex_collection::BiFMIndexCursor{index};
        REQUIRE(cursor.count() == index.size());
        REQUIRE(!cursor.empty());
        REQUIRE(cursor.lb == 0);
        REQUIRE(cursor.len == index.size());

        SECTION("case 1 - left") {
            auto c2 = cursor.extendLeft(1);
            REQUIRE(c2.len == 1);
            auto [e, steps] = index.locate(c2.lb);
            CHECK(steps == 0);
            auto [seq, pos] = e;
            CHECK(seq == 0);
            CHECK(pos == 0);
        }
        SECTION("case 2 - left") {
            auto c2 = cursor.extendLeft(2);
            REQUIRE(c2.len == 1);
            auto [e, steps] = index.locate(c2.lb);
            CHECK(steps == 0);
            auto [seq, pos] = e;
            CHECK(seq == 0);
            CHECK(pos == 1);
        }
        SECTION("case 0 - left") {
            auto c2 = cursor.extendLeft(0);
            REQUIRE(c2.len == 1);
            {
                auto [e, steps] = index.locate(c2.lb);
                CHECK(steps == 0);
                auto [seq, pos] = e;
                CHECK(seq == 0);
                CHECK(pos == 2);
            }
        }

        SECTION("case 1 - right") {
            auto c2 = cursor.extendRight(1);
            REQUIRE(c2.len == 1);
            auto [e, steps] = index.locate(c2.lb);
            CHECK(steps == 0);
            auto [seq, pos] = e;
            CHECK(seq == 0);
            CHECK(pos == 0);
        }
        SECTION("case 2 - right") {
            auto c2 = cursor.extendRight(2);
            REQUIRE(c2.len == 1);
            auto [e, steps] = index.locate(c2.lb);
            CHECK(steps == 0);
            auto [seq, pos] = e;
            CHECK(seq == 0);
            CHECK(pos == 1);
        }
        SECTION("case 0 - right") {
            auto c2 = cursor.extendRight(0);
            REQUIRE(c2.len == 1);
            {
                auto [e, steps] = index.locate(c2.lb);
                CHECK(steps == 0);
                auto [seq, pos] = e;
                CHECK(seq == 0);
                CHECK(pos == 2);
            }
        }
    }
}

TEST_CASE("checking bidirectional fm index without delimiters - bwt/sa", "[bifmindex-nd]") {
    auto input = std::vector<std::vector<uint8_t>> {std::vector<uint8_t>{1, 2, 0, 0, 1, 2}};
    auto expectedSA = std::vector<std::tuple<std::tuple<size_t, size_t>, size_t>> {
        {{0, 2}, 0},
        {{0, 3}, 0},
        {{0, 0}, 0},
        {{0, 4}, 0},
        {{0, 1}, 0},
        {{0, 5}, 0},
    };

    auto index = fmindex_collection::BiFMIndex<255>{input, /*.samplingRate=*/ 1, /*.threadNbr=*/1, false};
    CHECK(index.size() == expectedSA.size());

    for (size_t i{}; i < index.size(); ++i) {
        CHECK(index.locate(i) == expectedSA[i]);
    }
}
