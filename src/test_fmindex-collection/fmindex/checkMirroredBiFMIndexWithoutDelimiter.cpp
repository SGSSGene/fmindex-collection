// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "../string/allStrings.h"
#include "../string/utils.h"

#include <catch2/catch_all.hpp>
#include <fmindex-collection/fmindex/BiFMIndex.h>
#include <fmindex-collection/fmindex/BiFMIndexCursor.h>
#include <fstream>

TEST_CASE("checking mirrored bidirectional fm index without delimiters", "[mirroredbifmindex-nd]") {
    using Index = fmc::BiFMIndex<255>::NoDelim::ReuseRev;
    auto input = std::vector<std::vector<uint8_t>> {std::vector<uint8_t>{1, 2, 0}};
    SECTION("test index without delimiter") {
        auto index = Index{input, /*.samplingRate=*/ 1, /*.threadNbr=*/1, /*.seqOffset=*/0, /*.mirrorInput=*/true};
        auto cursor = fmc::BiFMIndexCursor{index};

        REQUIRE(index.size() == 6);
        REQUIRE(cursor.count() == index.size());
        REQUIRE(!cursor.empty());
        REQUIRE(cursor.lb == 0);
        REQUIRE(cursor.len == index.size());

        CHECK(index.bwt.symbol(0) == 2);
        CHECK(index.bwt.symbol(1) == 0);
        CHECK(index.bwt.symbol(2) == 2);
        CHECK(index.bwt.symbol(3) == 1);
        CHECK(index.bwt.symbol(4) == 1);
        CHECK(index.bwt.symbol(5) == 0);
        CHECK(std::get<1>(index.locate(0)) == 2);
        CHECK(std::get<1>(index.locate(1)) == 2);
        CHECK(std::get<1>(index.locate(2)) == 0);
        CHECK(std::get<1>(index.locate(3)) == 0);
        CHECK(std::get<1>(index.locate(4)) == 1);
        CHECK(std::get<1>(index.locate(5)) == 1);

        SECTION("case 1 - left") {
            auto c2 = cursor.extendLeft(1);
            REQUIRE(c2.len == 2);
            {
                auto [seq, pos, offset] = index.locate(c2.lb);
                CHECK(seq == 1);
                CHECK(pos+offset == 0);
            }
            {
                auto [seq, pos, offset] = index.locate(c2.lb+1);
                CHECK(seq == 0);
                CHECK(pos+offset == 0);
            }
        }
        SECTION("case 2 - left") {
            auto c2 = cursor.extendLeft(2);
            REQUIRE(c2.len == 2);
            {
                auto [seq, pos, offset] = index.locate(c2.lb);
                CHECK(seq == 0);
                CHECK(pos+offset == 1);
            }
            {
                auto [seq, pos, offset] = index.locate(c2.lb+1);
                CHECK(seq == 1);
                CHECK(pos+offset == 1);
            }
        }
        SECTION("case 0 - left") {
            auto c2 = cursor.extendLeft(0);
            REQUIRE(c2.len == 2);
            {
                auto [seq, pos, offset] = index.locate(c2.lb);
                CHECK(seq == 0);
                CHECK(pos+offset == 2);
            }
            {
                auto [seq, pos, offset] = index.locate(c2.lb+1);
                CHECK(seq == 1);
                CHECK(pos+offset == 2);
            }
        }

        SECTION("case 1 - right") {
            auto c2 = cursor.extendRight(1);
            REQUIRE(c2.len == 2);
            {
                auto [seq, pos, offset] = index.locate(c2.lb);
                CHECK(seq == 1);
                CHECK(pos+offset == 0);
            }
            {
                auto [seq, pos, offset] = index.locate(c2.lb+1);
                CHECK(seq == 0);
                CHECK(pos+offset == 0);
            }
        }
        SECTION("case 2 - right") {
            auto c2 = cursor.extendRight(2);
            REQUIRE(c2.len == 2);
            {
                auto [seq, pos, offset] = index.locate(c2.lb);
                CHECK(seq == 0);
                CHECK(pos+offset == 1);
            }
            {
                auto [seq, pos, offset] = index.locate(c2.lb+1);
                CHECK(seq == 1);
                CHECK(pos+offset == 1);
            }
        }
        SECTION("case 0 - right") {
            auto c2 = cursor.extendRight(0);
            REQUIRE(c2.len == 2);
            {
                auto [seq, pos, offset] = index.locate(c2.lb);
                CHECK(seq == 0);
                CHECK(pos+offset == 2);
            }
            {
                auto [seq, pos, offset] = index.locate(c2.lb+1);
                CHECK(seq == 1);
                CHECK(pos+offset == 2);
            }
        }
    }
}
