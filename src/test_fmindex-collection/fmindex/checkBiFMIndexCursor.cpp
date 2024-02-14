// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include "../occtables/allTables.h"

#include <catch2/catch_all.hpp>
#include <fmindex-collection/fmindex/BiFMIndex.h>
#include <fmindex-collection/fmindex/BiFMIndexCursor.h>

TEST_CASE("checking bidirectional fm index cursor", "[BiFMIndexCursor]") {

    auto data = std::vector<std::vector<uint8_t>>{std::vector<uint8_t>{1, 1, 1, 1, 2, 2, 2}};
    using OccTable = fmindex_collection::occtable::Bitvector<256>;
    using Index = fmindex_collection::BiFMIndex<OccTable>;
    auto index = Index{data, 1, 1};

    auto cursor = fmindex_collection::BiFMIndexCursor{index};
    REQUIRE(cursor.count() == index.size());
    REQUIRE(!cursor.empty());
    REQUIRE(cursor.lb == 0);
    REQUIRE(cursor.len == index.size());

    SECTION("extending to the left") {
        SECTION("check if cursor can be extended with 0") {
            auto cursor2 = cursor.extendLeft(0);
            REQUIRE(cursor2.count() == 1);
            REQUIRE(cursor2.lb == 0);
        }

        SECTION("check if cursor can be extended with 1") {
            auto cursor2 = cursor.extendLeft(1);
            REQUIRE(cursor2.count() == 4);
            REQUIRE(cursor2.lb == 1);
        }

        SECTION("check if cursor can be extended with 2") {
            auto cursor2 = cursor.extendLeft(2);
            REQUIRE(cursor2.count() == 3);
            REQUIRE(cursor2.lb == 5);
        }

        SECTION("check if cursor can be extended with 3") {
            auto cursor2 = cursor.extendLeft(3);
            REQUIRE(cursor2.count() == 0);
            REQUIRE(cursor2.lb == 8);
        }

        SECTION("check all symbols can be extended") {
            auto allCursor = cursor.extendLeft();
            for (auto i{0}; i < 256; ++i) {
                auto cursor2 = cursor.extendLeft(i);
                CHECK(cursor2.index == allCursor[i].index);
                CHECK(cursor2.lb == allCursor[i].lb);
                CHECK(cursor2.len == allCursor[i].len);
            }
        }
    }

    SECTION("extending to the right") {
        SECTION("check if cursor can be extended with 0") {
            auto cursor2 = cursor.extendRight(0);
            REQUIRE(cursor2.count() == 1);
            REQUIRE(cursor2.lb == 0);
        }

        SECTION("check if cursor can be extended with 1") {
            auto cursor2 = cursor.extendRight(1);
            REQUIRE(cursor2.count() == 4);
            REQUIRE(cursor2.lb == 1);
        }

        SECTION("check if cursor can be extended with 2") {
            auto cursor2 = cursor.extendRight(2);
            REQUIRE(cursor2.count() == 3);
            REQUIRE(cursor2.lb == 5);
        }

        SECTION("check if cursor can be extended with 3") {
            auto cursor2 = cursor.extendRight(3);
            REQUIRE(cursor2.count() == 0);
            REQUIRE(cursor2.lb == 8);
        }

        SECTION("check all symbols can be extended") {
            auto allCursor = cursor.extendRight();
            for (auto i{0}; i < 256; ++i) {
                auto cursor2 = cursor.extendRight(i);
                CHECK(cursor2.index == allCursor[i].index);
                CHECK(cursor2.lb == allCursor[i].lb);
                CHECK(cursor2.len == allCursor[i].len);
            }
        }
    }

    SECTION("check if iteration works") {
        auto entries = std::unordered_set<size_t>{0, 1, 2, 3, 4, 5, 6, 7};
        for (auto pos : cursor) {
            REQUIRE(entries.contains(pos));
            entries.erase(pos);
        }
        CHECK(entries.empty());
    }
}
