// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#include "allTables.h"

#include <fmindex-collection/BiFMIndex.h>
#include <fmindex-collection/search/SearchPseudo.h>
#include <search_schemes/generator/all.h>
#include <catch2/catch.hpp>

TEMPLATE_TEST_CASE("searching with PseudoSearch", "[search]", ALLTABLES) {
    using OccTable = TestType;

    auto input  = std::vector<uint8_t>{'A', 'A', 'A', 'C', 'A', 'A', 'A', 'C', 'A', 'A', 'A'};

    auto index = fmindex_collection::BiFMIndex<OccTable>{std::vector<std::vector<uint8_t>>{input}, /*samplingRate*/1, /*threadNbr*/1};

    SECTION("check symbol call to occurrence table") {
        REQUIRE(input.size()+1 == index.size());
        CHECK(index.occ.symbol( 0) == 'A');
        CHECK(index.occ.symbol( 1) == 'A');
        CHECK(index.occ.symbol( 2) == 'A');
        CHECK(index.occ.symbol( 3) == 'C');
        CHECK(index.occ.symbol( 4) == 'C');
        CHECK(index.occ.symbol( 5) == '\0');
        CHECK(index.occ.symbol( 6) == 'A');
        CHECK(index.occ.symbol( 7) == 'A');
        CHECK(index.occ.symbol( 8) == 'A');
        CHECK(index.occ.symbol( 9) == 'A');
        CHECK(index.occ.symbol(10) == 'A');
        CHECK(index.occ.symbol(11) == 'A');
    }

    auto query = std::vector<std::vector<uint8_t>> {std::vector<uint8_t>{'A'}};
    auto search_scheme = search_schemes::generator::backtracking(1, 0, 0);
    fmindex_collection::search_pseudo::search<true>(index, query, search_scheme, [](auto qidx, auto result, auto errors) {
        CHECK(qidx == 0);
        CHECK(errors == 0);
        CHECK(result.lb == 1);
        CHECK(result.count() == 9);
    });

}
TEMPLATE_TEST_CASE("searching with collection and PseudoSearch", "[collection]", ALLTABLES) {
    using OccTable = TestType;

    auto input  = std::vector<std::vector<uint8_t>>{{'A', 'A', 'A', 'C', 'A', 'A', 'A', 'C', 'A', 'A', 'A'},
                                                    {'A', 'A', 'A', 'B', 'A', 'A', 'A', 'B', 'A', 'A', 'A'}};

    auto index = fmindex_collection::BiFMIndex<OccTable>{input, /*samplingRate*/1, /*threadNbr*/1};

    SECTION("check symbol call to occurrence table") {
        auto expected = std::vector<uint8_t>{'A', 'A', 'A', 'A', 'A', 'A', 'B', 'C', 'B', '\0', 'C', '\0',
                                             'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'};
        REQUIRE(index.size() == expected.size());
        for (size_t i{0}; i < expected.size(); ++i) {
            INFO(i);
            CHECK(index.occ.symbol(i) == expected[i]);
        }
    }

    auto query = std::vector<std::vector<uint8_t>> {std::vector<uint8_t>{'A'}};
    auto search_scheme = search_schemes::generator::backtracking(1, 0, 0);
    fmindex_collection::search_pseudo::search<true>(index, query, search_scheme, [](auto qidx, auto result, auto errors) {
        CHECK(qidx == 0);
        CHECK(errors == 0);
        CHECK(result.lb == 2);
        CHECK(result.count() == 18);
    });

    auto expected = std::vector<std::tuple<size_t, size_t>> {
        std::make_tuple(1ull, 11ull),
        std::make_tuple(0ull, 11ull),
        std::make_tuple(1ull, 10ull),
        std::make_tuple(0ull, 10ull),
        std::make_tuple(1ull,  9ull),
        std::make_tuple(0ull,  9ull),
        std::make_tuple(1ull,  8ull),
        std::make_tuple(0ull,  8ull),
        std::make_tuple(1ull,  4ull),
        std::make_tuple(1ull,  0ull),
        std::make_tuple(0ull,  4ull),
        std::make_tuple(0ull,  0ull),
        std::make_tuple(1ull,  5ull),
        std::make_tuple(1ull,  1ull),
        std::make_tuple(0ull,  5ull),
        std::make_tuple(0ull,  1ull),
        std::make_tuple(1ull,  6ull),
        std::make_tuple(1ull,  2ull),
        std::make_tuple(0ull,  6ull),
        std::make_tuple(0ull,  2ull),
        std::make_tuple(1ull,  7ull),
        std::make_tuple(1ull,  3ull),
        std::make_tuple(0ull,  7ull),
        std::make_tuple(0ull,  3ull)
    };

    for (size_t i{0}; i < expected.size(); ++i) {
        INFO(i);
        auto [il, pl] = index.locate(i);
        auto [ir, pr] = expected[i];
        CHECK(il == ir);
        CHECK(pl == pr);
    }
}
