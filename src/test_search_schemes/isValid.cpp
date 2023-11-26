// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include <catch2/catch_all.hpp>
#include <search_schemes/isValid.h>

namespace ss = search_schemes;

TEST_CASE("check is valid", "[isValid]") {

    auto search = ss::Search{
        {0},
        {0},
        {0}
    };
    CHECK(ss::isValid(search));

    CHECK(ss::isValid(ss::Search {
        {0, 1},
        {0, 0},
        {0, 0},
    }));

    CHECK(ss::isValid(ss::Search {
        {1, 0},
        {0, 0},
        {0, 0},
    }));
    CHECK(ss::isValid(ss::Search {
        {0, 1, 2},
        {0, 0, 0},
        {0, 0, 0},
    }));
    CHECK(ss::isValid(ss::Search {
        {1, 0, 2},
        {0, 0, 0},
        {0, 0, 0},
    }));
    CHECK(ss::isValid(ss::Search {
        {1, 2, 0},
        {0, 0, 0},
        {0, 0, 0},
    }));
    CHECK(ss::isValid(ss::Search {
        {2, 1, 0},
        {0, 0, 0},
        {0, 0, 0},
    }));
    CHECK(not ss::isValid(ss::Search {
        {0, 2, 1},
        {0, 0, 0},
        {0, 0, 0},
    }));
    CHECK(not ss::isValid(ss::Search {
        {2, 0, 1},
        {0, 0, 0},
        {0, 0, 0},
    }));
    CHECK(not ss::isValid(ss::Search {
        {0, 0, 2},
        {0, 0, 0},
        {0, 0, 0},
    }));
}
