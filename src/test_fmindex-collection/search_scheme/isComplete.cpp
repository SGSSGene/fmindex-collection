// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <catch2/catch_all.hpp>
#include <fmindex-collection/search_scheme/isComplete.h>

namespace ss = fmc::search_scheme;

TEST_CASE("check is complete", "[isComplete]") {

    CHECK(ss::isComplete(ss::Scheme{ss::Search {
        {0, 1},
        {0, 0},
        {0, 0},
    }}, 0, 0));

    CHECK(not ss::isComplete(ss::Scheme{ss::Search {
        {0, 1},
        {0, 0},
        {0, 0},
    }}, 0, 1));

    CHECK(ss::isComplete(ss::Scheme{ss::Search {
        {0, 1},
        {0, 0},
        {1, 1},
    }}, 0, 1));

    CHECK(ss::isComplete(ss::Scheme{ss::Search {
        {0, 1},
        {0, 1},
        {1, 1},
    }}, 1, 1));
}
