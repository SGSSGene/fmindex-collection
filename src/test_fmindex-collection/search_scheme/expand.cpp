// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <catch2/catch_all.hpp>
#include <fmindex-collection/search_scheme/expand.h>
#include <fmindex-collection/search_scheme/isValid.h>

namespace ss = fmindex_collection::search_scheme;

TEST_CASE("check expand", "[expand]") {
    SECTION("no errors 2 to 10") {
        auto expected = ss::Scheme{ss::Search {
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        }};

        auto real = ss::expand(ss::Scheme{ss::Search {
            {0, 1},
            {0, 0},
            {0, 0},
        }}, 10);
        CHECK(ss::isValid(real));
        CHECK(expected == real);
    }

    SECTION("1 error 2 to 4") {
        auto expected = ss::Scheme{ss::Search {
            {0, 1, 2, 3},
            {0, 0, 0, 0},
            {0, 0, 1, 1},
        }};

        auto real = ss::expand(ss::Scheme{ss::Search {
            {0, 1},
            {0, 0},
            {0, 1},
        }}, 4);
        CHECK(ss::isValid(real));
        CHECK(expected == real);
    }
    SECTION("1 error 2 to 3") {
        auto expected = ss::Scheme{ss::Search {
            {0, 1, 2},
            {0, 0, 0},
            {0, 0, 1},
        }};

        auto real = ss::expand(ss::Scheme{ss::Search {
            {0, 1},
            {0, 0},
            {0, 1},
        }}, 3);
        CHECK(ss::isValid(real));
        CHECK(expected == real);
    }


}
