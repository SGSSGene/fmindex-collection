// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: CC0-1.0

#include "../string/allStrings.h"
#include "../string/utils.h"

#include <catch2/catch_all.hpp>
#include <fmindex-collection/fmindex/BiFMIndexKStep.h>
#include <fmindex-collection/fmindex/BiFMIndexKStepCursor.h>

TEST_CASE("checking bidirectional fm index with kstep cursor", "[bifmindexcursor_kstep]") {
    auto data = std::vector<std::vector<uint8_t>>{std::vector<uint8_t>{1, 1, 1, 1, 2, 2, 2}};
    using Index = fmc::BiFMIndexKStep<4>;
    auto index = Index{data, /*.samplingRate=*/1, /*.threadNbr=*/1};

    auto cursor = fmc::BiFMIndexKStepCursor{index};
    REQUIRE(cursor.count() == index.size());
    REQUIRE(!cursor.empty());
    REQUIRE(cursor.lb == 0);
    REQUIRE(cursor.len == index.size());

    SECTION("extending to the left") {
        SECTION("check if cursor can be extended with 0") {
            auto cursor2 = cursor.extendLeft(0);
            REQUIRE(cursor2.count() == 1);
            REQUIRE(cursor2.lb == 0);
            REQUIRE(cursor2.steps == 1);
        }

        SECTION("check if cursor can be extended with 1") {
            auto cursor2 = cursor.extendLeft(1);
            REQUIRE(cursor2.count() == 4);
            REQUIRE(cursor2.lb == 1);
            REQUIRE(cursor2.steps == 1);
        }

        SECTION("check if cursor can be extended with 2") {
            auto cursor2 = cursor.extendLeft(2);
            REQUIRE(cursor2.count() == 3);
            REQUIRE(cursor2.lb == 5);
            REQUIRE(cursor2.steps == 1);
        }

        SECTION("check if cursor can be extended with 3") {
            auto cursor2 = cursor.extendLeft(3);
            REQUIRE(cursor2.count() == 0);
            REQUIRE(cursor2.lb == 8);
            REQUIRE(cursor2.steps == 1);
        }

        SECTION("check all symbols can be extended") {
            auto allCursor = cursor.extendLeft();
            for (auto i{0}; i < 4; ++i) {
                auto cursor2 = cursor.extendLeft(i);
                CHECK(cursor2.index == allCursor[i].index);
                CHECK(cursor2.lb == allCursor[i].lb);
                CHECK(cursor2.len == allCursor[i].len);
                REQUIRE(cursor2.steps == allCursor[i].steps);
            }
        }

        SECTION("check if cursor can be extended twice, searching for 12") {
            auto cursor2 = cursor.extendLeft(2).extendLeft(1);
            CHECK(cursor2.count() == 1);
            CHECK(cursor2.lb == 4);
            CHECK(cursor2.lbRev == 5);
            CHECK(cursor2.steps == 2);
        }

        SECTION("check if cursor can be extended twice, searching for 12") {
            auto symbs = std::array<size_t, 2>{2, 1};
            auto cursor2 = cursor.extendLeftKStep(symbs);
            CHECK(cursor2.count() == 1);
            CHECK(cursor2.lb == 4);
            CHECK(cursor2.lbRev == 5);
            CHECK(cursor2.steps == 2);
        }
    }

    SECTION("extending to the right") {
        SECTION("check if cursor can be extended with 0") {
            auto cursor2 = cursor.extendRight(0);
            REQUIRE(cursor2.count() == 1);
            REQUIRE(cursor2.lb == 0);
            REQUIRE(cursor2.steps == 1);
        }

        SECTION("check if cursor can be extended with 1") {
            auto cursor2 = cursor.extendRight(1);
            REQUIRE(cursor2.count() == 4);
            REQUIRE(cursor2.lb == 1);
            REQUIRE(cursor2.steps == 1);
        }

        SECTION("check if cursor can be extended with 2") {
            auto cursor2 = cursor.extendRight(2);
            REQUIRE(cursor2.count() == 3);
            REQUIRE(cursor2.lb == 5);
            REQUIRE(cursor2.steps == 1);
        }

        SECTION("check if cursor can be extended with 3") {
            auto cursor2 = cursor.extendRight(3);
            REQUIRE(cursor2.count() == 0);
            REQUIRE(cursor2.lb == 8);
            REQUIRE(cursor2.steps == 1);
        }

        SECTION("check all symbols can be extended") {
            auto allCursor = cursor.extendRight();
            for (auto i{0}; i < 4; ++i) {
                auto cursor2 = cursor.extendRight(i);
                CHECK(cursor2.index == allCursor[i].index);
                CHECK(cursor2.lb == allCursor[i].lb);
                CHECK(cursor2.len == allCursor[i].len);
                REQUIRE(cursor2.steps == allCursor[i].steps);
            }
        }

        SECTION("check if cursor can be extended twice, searching for 21") {
            auto cursor2 = cursor.extendRight(1).extendRight(2);
            CHECK(cursor2.count() == 1);
            CHECK(cursor2.lb == 4);
            CHECK(cursor2.lbRev == 5);
            CHECK(cursor2.steps == 2);
        }

        SECTION("check if cursor can be extended twice, searching for 21") {
            auto symbs = std::array<size_t, 2>{1, 2};
            auto cursor2 = cursor.extendRightKStep(symbs);
            CHECK(cursor2.count() == 1);
            CHECK(cursor2.lb == 4);
            CHECK(cursor2.lbRev == 5);
            CHECK(cursor2.steps == 2);
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
