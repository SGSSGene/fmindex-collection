// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "../string/allStrings.h"

#include <catch2/catch_all.hpp>
#include <fmindex-collection/locate.h>
#include <fmindex-collection/search/Backtracking.h>
#include <fmindex-collection/search_scheme/generator/all.h>

TEST_CASE("searching with backtracking", "[search]") {
    auto input  = std::vector<uint8_t>{'A', 'A', 'A', 'C', 'A', 'A', 'A', 'C', 'A', 'A', 'A'};

    auto index = fmc::BiFMIndex<255>{std::vector<std::vector<uint8_t>>{input}, /*samplingRate*/1, /*threadNbr*/1};

    SECTION("check symbol call to occurrence table") {
        REQUIRE(input.size()+1 == index.size());
        CHECK(index.bwt.symbol( 0) == 'A');
        CHECK(index.bwt.symbol( 1) == 'A');
        CHECK(index.bwt.symbol( 2) == 'A');
        CHECK(index.bwt.symbol( 3) == 'C');
        CHECK(index.bwt.symbol( 4) == 'C');
        CHECK(index.bwt.symbol( 5) == '\0');
        CHECK(index.bwt.symbol( 6) == 'A');
        CHECK(index.bwt.symbol( 7) == 'A');
        CHECK(index.bwt.symbol( 8) == 'A');
        CHECK(index.bwt.symbol( 9) == 'A');
        CHECK(index.bwt.symbol(10) == 'A');
        CHECK(index.bwt.symbol(11) == 'A');
    }

    auto query = std::vector<std::vector<uint8_t>> {std::vector<uint8_t>{'A'}};
    fmc::search_backtracking::search(index, query, 0, [](auto qidx, auto result, auto errors) {
        CHECK(qidx == 0);
        CHECK(errors == 0);
        CHECK(result.lb == 1);
        CHECK(result.count() == 9);
    });
}

TEST_CASE("searching with collection and backtracking", "[collection]") {
    auto input  = std::vector<std::vector<uint8_t>>{{'A', 'A', 'A', 'C', 'A', 'A', 'A', 'C', 'A', 'A', 'A'},
                                                    {'A', 'A', 'A', 'B', 'A', 'A', 'A', 'B', 'A', 'A', 'A'}};

    auto index = fmc::BiFMIndex<255>{input, /*samplingRate*/1, /*threadNbr*/1};

    SECTION("check symbol call to occurrence table") {
        auto expected = std::vector<uint8_t>{'A', 'A', 'A', 'A', 'A', 'A', 'B', 'C', 'B', '\0', 'C', '\0',
                                             'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'};
        REQUIRE(index.size() == expected.size());
        for (size_t i{0}; i < expected.size(); ++i) {
            INFO(i);
            CHECK(index.bwt.symbol(i) == expected[i]);
        }
    }

    auto query = std::vector<std::vector<uint8_t>> {std::vector<uint8_t>{'A'}};
    fmc::search_backtracking::search(index, query, 0, [](auto qidx, auto result, auto errors) {
        CHECK(qidx == 0);
        CHECK(errors == 0);
        CHECK(result.lb == 2);
        CHECK(result.count() == 18);
    });

    auto expected = std::vector<std::tuple<uint32_t, uint32_t>> {
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
        auto [entry, offset] = index.locate(i);
        auto [il, pl] = entry;
        pl += offset;

        auto [ir, pr] = expected[i];
        CHECK(il == ir);
        CHECK(pl == pr);
    }
}

TEST_CASE("searching with backtracking with FMIndex", "[search]") {
    auto input  = std::vector<uint8_t>{'A', 'A', 'A', 'C', 'A', 'A', 'A', 'C', 'A', 'A', 'A'};

    auto index = fmc::FMIndex<255>{input, /*samplingRate*/1, /*samplingRate*/1};

    SECTION("check symbol call to occurrence table") {
        REQUIRE(input.size()+1 == index.size());
        CHECK(index.bwt.symbol( 0) == 'A');
        CHECK(index.bwt.symbol( 1) == 'A');
        CHECK(index.bwt.symbol( 2) == 'A');
        CHECK(index.bwt.symbol( 3) == 'C');
        CHECK(index.bwt.symbol( 4) == 'C');
        CHECK(index.bwt.symbol( 5) == '\0');
        CHECK(index.bwt.symbol( 6) == 'A');
        CHECK(index.bwt.symbol( 7) == 'A');
        CHECK(index.bwt.symbol( 8) == 'A');
        CHECK(index.bwt.symbol( 9) == 'A');
        CHECK(index.bwt.symbol(10) == 'A');
        CHECK(index.bwt.symbol(11) == 'A');
    }

    auto query = std::vector<std::vector<uint8_t>> {std::vector<uint8_t>{'A'}};
    fmc::search_backtracking::search(index, query, 0, [](auto qidx, auto result, auto errors) {
        CHECK(qidx == 0);
        CHECK(errors == 0);
        CHECK(result.lb == 1);
        CHECK(result.count() == 9);
    });
}

TEST_CASE("searching with collection and backtracking with FMIndex", "[collection]") {
    auto input  = std::vector<std::vector<uint8_t>>{{'A', 'A', 'A', 'C', 'A', 'A', 'A', 'C', 'A', 'A', 'A'},
                                                    {'A', 'A', 'A', 'B', 'A', 'A', 'A', 'B', 'A', 'A', 'A'}};

    auto index = fmc::FMIndex<255>{input, /*samplingRate*/1, /*threadNbr*/1};

    SECTION("check symbol call to occurrence table") {
        auto expected = std::vector<uint8_t>{'A', 'A', 'A', 'A', 'A', 'A', 'B', 'C', 'B', '\0', 'C', '\0',
                                             'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'};
        REQUIRE(index.size() == expected.size());
        for (size_t i{0}; i < expected.size(); ++i) {
            INFO(i);
            CHECK(index.bwt.symbol(i) == expected[i]);
        }
    }

    auto query = std::vector<std::vector<uint8_t>> {std::vector<uint8_t>{'A'}};
    fmc::search_backtracking::search(index, query, 0, [](auto qidx, auto result, auto errors) {
        CHECK(qidx == 0);
        CHECK(errors == 0);
        CHECK(result.lb == 2);
        CHECK(result.count() == 18);
    });

    auto expected = std::vector<std::tuple<uint32_t, uint32_t>> {
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
        auto [entry, offset] = index.locate(i);
        auto [il, pl] = entry;
        pl += offset;

        auto [ir, pr] = expected[i];
        CHECK(il == ir);
        CHECK(pl == pr);
    }
}

TEST_CASE("searching with backtracking with ReverseFMIndex", "[search]") {
    using String = fmc::string::InterleavedBitvector16<255>;

    auto input  = std::vector<uint8_t>{'A', 'A', 'A', 'C', 'A', 'A', 'A', 'C', 'A', 'A', 'A'};

    auto index = fmc::ReverseFMIndex<String>{std::vector<std::vector<uint8_t>>{input}, /*samplingRate*/1, /*threadNbr*/1};

    SECTION("check symbol call to occurrence table") {
        REQUIRE(input.size()+1 == index.size());
        CHECK(index.bwt.symbol( 0) == 'A');
        CHECK(index.bwt.symbol( 1) == 'A');
        CHECK(index.bwt.symbol( 2) == 'A');
        CHECK(index.bwt.symbol( 3) == 'C');
        CHECK(index.bwt.symbol( 4) == 'C');
        CHECK(index.bwt.symbol( 5) == '\0');
        CHECK(index.bwt.symbol( 6) == 'A');
        CHECK(index.bwt.symbol( 7) == 'A');
        CHECK(index.bwt.symbol( 8) == 'A');
        CHECK(index.bwt.symbol( 9) == 'A');
        CHECK(index.bwt.symbol(10) == 'A');
        CHECK(index.bwt.symbol(11) == 'A');
    }

    auto query = std::vector<std::vector<uint8_t>> {std::vector<uint8_t>{'A'}};
    fmc::search_backtracking::search(index, query, 0, [](auto qidx, auto result, auto errors) {
        CHECK(qidx == 0);
        CHECK(errors == 0);
        CHECK(result.lb == 1);
        CHECK(result.count() == 9);
    });

}
TEST_CASE("searching with collection and backtracking with ReverseFMIndex", "[collection]") {
    using String = fmc::string::InterleavedBitvector16<255>;

    auto input  = std::vector<std::vector<uint8_t>>{{'A', 'A', 'A', 'C', 'A', 'A', 'A', 'C', 'A', 'A', 'A'},
                                                    {'A', 'A', 'A', 'B', 'A', 'A', 'A', 'B', 'A', 'A', 'A'}};

    auto index = fmc::ReverseFMIndex<String>{input, /*samplingRate*/1, /*threadNbr*/1};

    SECTION("check symbol call to occurrence table") {
        auto expected = std::vector<uint8_t>{'A', 'A', 'A', 'A', 'A', 'A', 'B', 'C', 'B', '\0', 'C', '\0',
                                             'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'};
        REQUIRE(index.size() == expected.size());
        for (size_t i{0}; i < expected.size(); ++i) {
            INFO(i);
            CHECK(index.bwt.symbol(i) == expected[i]);
        }
    }

    auto query = std::vector<std::vector<uint8_t>> {std::vector<uint8_t>{'A'}};
    fmc::search_backtracking::search(index, query, 0, [](auto qidx, auto result, auto errors) {
        CHECK(qidx == 0);
        CHECK(errors == 0);
        CHECK(result.lb == 2);
        CHECK(result.count() == 18);
    });

    auto expected = std::vector<std::tuple<uint32_t, uint32_t>> {
        std::make_tuple(1ull,  0ull),
        std::make_tuple(0ull,  0ull),
        std::make_tuple(1ull,  1ull),
        std::make_tuple(0ull,  1ull),
        std::make_tuple(1ull,  2ull),
        std::make_tuple(0ull,  2ull),
        std::make_tuple(1ull,  3ull),
        std::make_tuple(0ull,  3ull),
        std::make_tuple(1ull,  7ull),
        std::make_tuple(1ull, 11ull),
        std::make_tuple(0ull,  7ull),
        std::make_tuple(0ull, 11ull),
        std::make_tuple(1ull,  6ull),
        std::make_tuple(1ull, 10ull),
        std::make_tuple(0ull,  6ull),
        std::make_tuple(0ull, 10ull),
        std::make_tuple(1ull,  5ull),
        std::make_tuple(1ull,  9ull),
        std::make_tuple(0ull,  5ull),
        std::make_tuple(0ull,  9ull),
        std::make_tuple(1ull,  4ull),
        std::make_tuple(1ull,  8ull),
        std::make_tuple(0ull,  4ull),
        std::make_tuple(0ull,  8ull)
    };

    for (size_t i{0}; i < expected.size(); ++i) {
        INFO(i);
        if (index.bwt.symbol(i) != 0) {
            auto [entry, offset] = index.locate(i);
            auto [il, pl] = entry;
            pl += offset;

            auto [ir, pr] = expected[i];
            CHECK(il == ir);
            CHECK(pl == pr);
        }
    }
}


TEST_CASE("backtracking with errors", "[collection]") {
    using Index = fmc::FMIndex<256>;

    auto input  = std::vector<std::vector<uint8_t>>{{'A', 'A', 'A', 'C', 'A', 'A', 'A', 'B', 'A', 'A', 'A'},
                                                    {'A', 'A', 'A', 'B', 'A', 'A', 'A', 'C', 'A', 'A', 'A'}};

    auto index = Index{input, /*samplingRate*/1, /*threadNbr*/1};

    auto queries = std::vector<std::vector<uint8_t>> {std::vector<uint8_t>{'C', 'C'}, std::vector<uint8_t>{'B', 'B'}};

    auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
    fmc::search_backtracking::search(index, queries, 1, [&](auto qidx, auto cursor, auto errors) {
        (void)errors;
        for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
            auto [sid, spos] = entry;
            results.emplace_back(qidx, sid, spos+offset);
        }
    });

    std::ranges::sort(results);

    auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
        {0, 0, 2},
        {0, 0, 3},
        {0, 1, 6},
        {0, 1, 7},
        {1, 0, 6},
        {1, 0, 7},
        {1, 1, 2},
        {1, 1, 3},
    };

    CHECK(results == expected);
}
