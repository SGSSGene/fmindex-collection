// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "../string/allStrings.h"

#include <catch2/catch_all.hpp>
#include <fmindex-collection/fmindex/FMIndex.h>
#include <fmindex-collection/fmindex/FMIndexCursor.h>
#include <fmindex-collection/suffixarray/CSA.h>
#include <fmindex-collection/suffixarray/utils.h>

TEST_CASE("checking unidirectional fm index cursor", "[fmindexcursor]") {

    auto data = std::vector<std::vector<uint8_t>>{std::vector<uint8_t>{1, 1, 1, 1, 2, 2, 2}};
    using Index = fmc::FMIndex<256>;
    auto index = Index{data, /*.samplingRate=*/1, /*.threadNbr=*/1};

    auto cursor = fmc::FMIndexCursor{index};
    REQUIRE(cursor.count() == index.size());
    REQUIRE(!cursor.empty());
    REQUIRE(cursor.lb == 0);
    REQUIRE(cursor.len == index.size());

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
    SECTION("check if iteration works") {
        auto entries = std::unordered_set<size_t>{0, 1, 2, 3, 4, 5, 6, 7};
        for (auto pos : cursor) {
            REQUIRE(entries.contains(pos));
            entries.erase(pos);
        }
        CHECK(entries.empty());
    }
}
