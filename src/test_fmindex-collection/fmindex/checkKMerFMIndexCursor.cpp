// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "../string/allStrings.h"

#include <catch2/catch_all.hpp>
#include <fmindex-collection/fmindex/KMerFMIndex.h>
#include <fmindex-collection/fmindex/KMerFMIndexCursor.h>

TEST_CASE("checking unidirectional kmer fm index cursor", "[kmerfmindexcursor]") {

    auto data = std::vector<uint8_t>{1, 1, 1, 1, 2, 2, 2};
    using String = fmc::string::MultiBitvector<256>;
    using Index = fmc::KMerFMIndex<String, 2>;
    auto index = Index{data, 1, 1};

    auto cursor = fmc::KMerFMIndexCursor{index};
    REQUIRE(cursor.count() == index.size());
    REQUIRE(!cursor.empty());
    REQUIRE(cursor.lb == 0);
    REQUIRE(cursor.len == index.size());
    REQUIRE(index.size() == 8);
    CHECK(index.bwt.symbol(0) == 2); // 0 111122 2
    CHECK(index.bwt.symbol(1) == 0); // 1 111222 0
    CHECK(index.bwt.symbol(2) == 1); // 1 112220 1
    CHECK(index.bwt.symbol(3) == 1); // 1 122201 1
    CHECK(index.bwt.symbol(4) == 1); // 1 222011 1
    CHECK(index.bwt.symbol(5) == 2); // 2 011112 2
    CHECK(index.bwt.symbol(6) == 2); // 2 201111 2
    CHECK(index.bwt.symbol(7) == 1); // 2 220111 1

    REQUIRE(index.kmerStarts.size() == 9);
    CHECK(index.kmerStarts.symbol(0) == 1);    // 0 111122 2    1
    CHECK(index.kmerStarts.symbol(1) == 1);    // 1 111222 0    1
    CHECK(index.kmerStarts.symbol(2) == 0);    // 1 112220 1    0
    CHECK(index.kmerStarts.symbol(3) == 0);    // 1 122201 1    0
    CHECK(index.kmerStarts.symbol(4) == 1);    // 1 222011 1    1
    CHECK(index.kmerStarts.symbol(5) == 1);    // 2 011112 2    1
    CHECK(index.kmerStarts.symbol(6) == 1);    // 2 201111 2    1
    CHECK(index.kmerStarts.symbol(7) == 0);    // 2 220111 1    0
    CHECK(index.kmerStarts.symbol(8) == 1);    //               1



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

    SECTION("check if we can continue cycling through the kmers") {
        auto c = cursor.clipToKMer();
        CHECK(c.lb == cursor.lb);
        CHECK(c.len == cursor.len);
        CHECK(c.index == cursor.index);

        // 0 111122 2    1
        // 1 111222 0    1
        // 1 112220 1    0
        // 1 122201 1    0
        // 1 222011 1    1
        // 2 011112 2    1
        // 2 201111 2    1
        // 2 220111 1    0
        //               1

        // pattern will be extended to length 1
        c = c.extendLeft(1).clipToKMer();
        CHECK(c.count() == 4);
        CHECK(c.lb == 1);

        // pattern will be extended to length 2
        c = c.extendLeft(1).clipToKMer();
        CHECK(c.count() == 3);
        CHECK(c.lb == 1);

        // pattern will be extended to length 3 and reduced back to length 2
        c = c.extendLeft(1).clipToKMer();
        CHECK(c.count() == 3);
        CHECK(c.lb == 1);
    }
}
