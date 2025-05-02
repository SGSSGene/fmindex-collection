// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include "../string/allStrings.h"

#include <catch2/catch_all.hpp>
#include <fmindex-collection/fmindex/MirroredBiFMIndex.h>
#include <fmindex-collection/fmindex/MirroredBiFMIndexCursor.h>
#include <fstream>

TEMPLATE_TEST_CASE("checking mirrored bidirectional fm index without delimiters", "[MirroredBiFMIndex-nd]", ALLRANKVECTORS(255)) {
    using OccTable = TestType;

    auto input = std::vector<std::vector<uint8_t>> {std::vector<uint8_t>{1, 2, 0}};
    SECTION("test index without delimiter") {
        auto index = fmindex_collection::MirroredBiFMIndex<OccTable>{input, /*.samplingRate=*/ 1, /*.threadNbr=*/1, false};

        auto cursor = fmindex_collection::MirroredBiFMIndexCursor{index};
        REQUIRE(cursor.count() == index.size());
        REQUIRE(!cursor.empty());
        REQUIRE(cursor.lb == 0);
        REQUIRE(cursor.len == index.size());

        SECTION("case 1 - left") {
            auto c2 = cursor.extendLeft(1);
            REQUIRE(c2.len == 2);
            {
                auto [seq, pos] = index.locate(c2.lb);
                CHECK(seq == 2);
                CHECK(pos == 2);
            }
            {
                auto [seq, pos] = index.locate(c2.lb+1);
                CHECK(seq == 1);
                CHECK(pos == 0);
            }
        }
        SECTION("case 2 - left") {
            auto c2 = cursor.extendLeft(2);
            REQUIRE(c2.len == 2);
            {
                auto [seq, pos] = index.locate(c2.lb);
                CHECK(seq == 1);
                CHECK(pos == 1);
            }
            {
                auto [seq, pos] = index.locate(c2.lb+1);
                CHECK(seq == 2);
                CHECK(pos == 1);
            }
        }
        SECTION("case 0 - left") {
            auto c2 = cursor.extendLeft(0);
            REQUIRE(c2.len == 6);
            {
                auto [seq, pos] = index.locate(c2.lb);
                CHECK(seq == 3);
                CHECK(pos == 1);
            }
            {
                auto [seq, pos] = index.locate(c2.lb+1);
                CHECK(seq == 3);
                CHECK(pos == 0);
            }
            {
                auto [seq, pos] = index.locate(c2.lb+2);
                CHECK(seq == 0);
                CHECK(pos == 0);
            }
            {
                auto [seq, pos] = index.locate(c2.lb+3);
                CHECK(seq == 1);
                CHECK(pos == 2);
            }
            {
                auto [seq, pos] = index.locate(c2.lb+4);
                CHECK(seq == 0);
                CHECK(pos == 1);
            }
            {
                auto [seq, pos] = index.locate(c2.lb+5);
                CHECK(seq == 2);
                CHECK(pos == 0);
            }
        }

        SECTION("case 1 - right") {
            auto c2 = cursor.extendRight(1);
            REQUIRE(c2.len == 2);
            {
                auto [seq, pos] = index.locate(c2.lb);
                CHECK(seq == 2);
                CHECK(pos == 2);
            }
            {
                auto [seq, pos] = index.locate(c2.lb+1);
                CHECK(seq == 1);
                CHECK(pos == 0);
            }
        }
        SECTION("case 2 - right") {
            auto c2 = cursor.extendRight(2);
            REQUIRE(c2.len == 2);
            {
                auto [seq, pos] = index.locate(c2.lb);
                CHECK(seq == 1);
                CHECK(pos == 1);
            }
            {
                auto [seq, pos] = index.locate(c2.lb+1);
                CHECK(seq == 2);
                CHECK(pos == 1);
            }
        }
        SECTION("case 0 - right") {
            auto c2 = cursor.extendRight(0);
            REQUIRE(c2.len == 6);
            {
                auto [seq, pos] = index.locate(c2.lb);
                CHECK(seq == 3);
                CHECK(pos == 1);
            }
            {
                auto [seq, pos] = index.locate(c2.lb+1);
                CHECK(seq == 3);
                CHECK(pos == 0);
            }
            {
                auto [seq, pos] = index.locate(c2.lb+2);
                CHECK(seq == 0);
                CHECK(pos == 0);
            }
            {
                auto [seq, pos] = index.locate(c2.lb+3);
                CHECK(seq == 1);
                CHECK(pos == 2);
            }
            {
                auto [seq, pos] = index.locate(c2.lb+4);
                CHECK(seq == 0);
                CHECK(pos == 1);
            }
            {
                auto [seq, pos] = index.locate(c2.lb+5);
                CHECK(seq == 2);
                CHECK(pos == 0);
            }
        }

    }
}
