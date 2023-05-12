// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#include <search_schemes/isComplete.h>
#include <catch2/catch.hpp>

namespace ss = search_schemes;

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
