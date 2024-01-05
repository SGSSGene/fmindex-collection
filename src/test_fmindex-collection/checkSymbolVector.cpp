// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include <catch2/catch_all.hpp>
#include <fmindex-collection/utils.h>

#include "allSymbolVectors.h"

TEMPLATE_TEST_CASE("check if rank on the symbol vectors is working", "[SymbolVector]", ALLSYMBOLVECTORS) {
    using Vector = TestType;
    INFO(typeid(Vector).name());

    auto text = std::vector<uint8_t>{'H', 'a', 'l', 'l', 'o', ' ', 'W', 'e', 'l', 't'};

    auto vec = Vector{std::span{text}};
    REQUIRE(vec.size() == text.size());

    SECTION("check that symbol() call works") {
        for (size_t i{0}; i < text.size(); ++i) {
            INFO(i);
            CHECK(vec.symbol(i) == text.at(i));
        }
    }
    CHECK(Vector::Sigma >= 128);


    SECTION("test complete vector on rank()") {
        CHECK(vec.rank( 0, ' ') == 0);
        CHECK(vec.rank( 1, ' ') == 0);
        CHECK(vec.rank( 2, ' ') == 0);
        CHECK(vec.rank( 3, ' ') == 0);
        CHECK(vec.rank( 4, ' ') == 0);
        CHECK(vec.rank( 5, ' ') == 0);
        CHECK(vec.rank( 6, ' ') == 1);
        CHECK(vec.rank( 7, ' ') == 1);
        CHECK(vec.rank( 8, ' ') == 1);
        CHECK(vec.rank( 9, ' ') == 1);
        CHECK(vec.rank(10, ' ') == 1);

        CHECK(vec.rank( 0, 'H') == 0);
        CHECK(vec.rank( 1, 'H') == 1);
        CHECK(vec.rank( 2, 'H') == 1);
        CHECK(vec.rank( 3, 'H') == 1);
        CHECK(vec.rank( 4, 'H') == 1);
        CHECK(vec.rank( 5, 'H') == 1);
        CHECK(vec.rank( 6, 'H') == 1);
        CHECK(vec.rank( 7, 'H') == 1);
        CHECK(vec.rank( 8, 'H') == 1);
        CHECK(vec.rank( 9, 'H') == 1);
        CHECK(vec.rank(10, 'H') == 1);

        CHECK(vec.rank( 0, 'W') == 0);
        CHECK(vec.rank( 1, 'W') == 0);
        CHECK(vec.rank( 2, 'W') == 0);
        CHECK(vec.rank( 3, 'W') == 0);
        CHECK(vec.rank( 4, 'W') == 0);
        CHECK(vec.rank( 5, 'W') == 0);
        CHECK(vec.rank( 6, 'W') == 0);
        CHECK(vec.rank( 7, 'W') == 1);
        CHECK(vec.rank( 8, 'W') == 1);
        CHECK(vec.rank( 9, 'W') == 1);
        CHECK(vec.rank(10, 'W') == 1);

        CHECK(vec.rank( 0, 'a') == 0);
        CHECK(vec.rank( 1, 'a') == 0);
        CHECK(vec.rank( 2, 'a') == 1);
        CHECK(vec.rank( 3, 'a') == 1);
        CHECK(vec.rank( 4, 'a') == 1);
        CHECK(vec.rank( 5, 'a') == 1);
        CHECK(vec.rank( 6, 'a') == 1);
        CHECK(vec.rank( 7, 'a') == 1);
        CHECK(vec.rank( 8, 'a') == 1);
        CHECK(vec.rank( 9, 'a') == 1);
        CHECK(vec.rank(10, 'a') == 1);

        CHECK(vec.rank( 0, 'e') == 0);
        CHECK(vec.rank( 1, 'e') == 0);
        CHECK(vec.rank( 2, 'e') == 0);
        CHECK(vec.rank( 3, 'e') == 0);
        CHECK(vec.rank( 4, 'e') == 0);
        CHECK(vec.rank( 5, 'e') == 0);
        CHECK(vec.rank( 6, 'e') == 0);
        CHECK(vec.rank( 7, 'e') == 0);
        CHECK(vec.rank( 8, 'e') == 1);
        CHECK(vec.rank( 9, 'e') == 1);
        CHECK(vec.rank(10, 'e') == 1);

        CHECK(vec.rank( 0, 'l') == 0);
        CHECK(vec.rank( 1, 'l') == 0);
        CHECK(vec.rank( 2, 'l') == 0);
        CHECK(vec.rank( 3, 'l') == 1);
        CHECK(vec.rank( 4, 'l') == 2);
        CHECK(vec.rank( 5, 'l') == 2);
        CHECK(vec.rank( 6, 'l') == 2);
        CHECK(vec.rank( 7, 'l') == 2);
        CHECK(vec.rank( 8, 'l') == 2);
        CHECK(vec.rank( 9, 'l') == 3);
        CHECK(vec.rank(10, 'l') == 3);

        CHECK(vec.rank( 0, 'o') == 0);
        CHECK(vec.rank( 1, 'o') == 0);
        CHECK(vec.rank( 2, 'o') == 0);
        CHECK(vec.rank( 3, 'o') == 0);
        CHECK(vec.rank( 4, 'o') == 0);
        CHECK(vec.rank( 5, 'o') == 1);
        CHECK(vec.rank( 6, 'o') == 1);
        CHECK(vec.rank( 7, 'o') == 1);
        CHECK(vec.rank( 8, 'o') == 1);
        CHECK(vec.rank( 9, 'o') == 1);
        CHECK(vec.rank(10, 'o') == 1);

        CHECK(vec.rank( 0, 't') == 0);
        CHECK(vec.rank( 1, 't') == 0);
        CHECK(vec.rank( 2, 't') == 0);
        CHECK(vec.rank( 3, 't') == 0);
        CHECK(vec.rank( 4, 't') == 0);
        CHECK(vec.rank( 5, 't') == 0);
        CHECK(vec.rank( 6, 't') == 0);
        CHECK(vec.rank( 7, 't') == 0);
        CHECK(vec.rank( 8, 't') == 0);
        CHECK(vec.rank( 9, 't') == 0);
        CHECK(vec.rank(10, 't') == 1);
    }

    SECTION("test complete vec 'H' for prefix_rank()") {
        CHECK(vec.prefix_rank( 0, ' ') == 0);
        CHECK(vec.prefix_rank( 1, ' ') == 0);
        CHECK(vec.prefix_rank( 2, ' ') == 0);
        CHECK(vec.prefix_rank( 3, ' ') == 0);
        CHECK(vec.prefix_rank( 4, ' ') == 0);
        CHECK(vec.prefix_rank( 5, ' ') == 0);
        CHECK(vec.prefix_rank( 6, ' ') == 1);
        CHECK(vec.prefix_rank( 7, ' ') == 1);
        CHECK(vec.prefix_rank( 8, ' ') == 1);
        CHECK(vec.prefix_rank( 9, ' ') == 1);
        CHECK(vec.prefix_rank(10, ' ') == 1);

        CHECK(vec.prefix_rank( 0, 'H') == 0);
        CHECK(vec.prefix_rank( 1, 'H') == 1);
        CHECK(vec.prefix_rank( 2, 'H') == 1);
        CHECK(vec.prefix_rank( 3, 'H') == 1);
        CHECK(vec.prefix_rank( 4, 'H') == 1);
        CHECK(vec.prefix_rank( 5, 'H') == 1);
        CHECK(vec.prefix_rank( 6, 'H') == 2);
        CHECK(vec.prefix_rank( 7, 'H') == 2);
        CHECK(vec.prefix_rank( 8, 'H') == 2);
        CHECK(vec.prefix_rank( 9, 'H') == 2);
        CHECK(vec.prefix_rank(10, 'H') == 2);

        CHECK(vec.prefix_rank( 0, 'W') == 0);
        CHECK(vec.prefix_rank( 1, 'W') == 1);
        CHECK(vec.prefix_rank( 2, 'W') == 1);
        CHECK(vec.prefix_rank( 3, 'W') == 1);
        CHECK(vec.prefix_rank( 4, 'W') == 1);
        CHECK(vec.prefix_rank( 5, 'W') == 1);
        CHECK(vec.prefix_rank( 6, 'W') == 2);
        CHECK(vec.prefix_rank( 7, 'W') == 3);
        CHECK(vec.prefix_rank( 8, 'W') == 3);
        CHECK(vec.prefix_rank( 9, 'W') == 3);
        CHECK(vec.prefix_rank(10, 'W') == 3);

        CHECK(vec.prefix_rank( 0, 'a') == 0);
        CHECK(vec.prefix_rank( 1, 'a') == 1);
        CHECK(vec.prefix_rank( 2, 'a') == 2);
        CHECK(vec.prefix_rank( 3, 'a') == 2);
        CHECK(vec.prefix_rank( 4, 'a') == 2);
        CHECK(vec.prefix_rank( 5, 'a') == 2);
        CHECK(vec.prefix_rank( 6, 'a') == 3);
        CHECK(vec.prefix_rank( 7, 'a') == 4);
        CHECK(vec.prefix_rank( 8, 'a') == 4);
        CHECK(vec.prefix_rank( 9, 'a') == 4);
        CHECK(vec.prefix_rank(10, 'a') == 4);

        CHECK(vec.prefix_rank( 0, 'e') == 0);
        CHECK(vec.prefix_rank( 1, 'e') == 1);
        CHECK(vec.prefix_rank( 2, 'e') == 2);
        CHECK(vec.prefix_rank( 3, 'e') == 2);
        CHECK(vec.prefix_rank( 4, 'e') == 2);
        CHECK(vec.prefix_rank( 5, 'e') == 2);
        CHECK(vec.prefix_rank( 6, 'e') == 3);
        CHECK(vec.prefix_rank( 7, 'e') == 4);
        CHECK(vec.prefix_rank( 8, 'e') == 5);
        CHECK(vec.prefix_rank( 9, 'e') == 5);
        CHECK(vec.prefix_rank(10, 'e') == 5);

        CHECK(vec.prefix_rank( 0, 'l') == 0);
        CHECK(vec.prefix_rank( 1, 'l') == 1);
        CHECK(vec.prefix_rank( 2, 'l') == 2);
        CHECK(vec.prefix_rank( 3, 'l') == 3);
        CHECK(vec.prefix_rank( 4, 'l') == 4);
        CHECK(vec.prefix_rank( 5, 'l') == 4);
        CHECK(vec.prefix_rank( 6, 'l') == 5);
        CHECK(vec.prefix_rank( 7, 'l') == 6);
        CHECK(vec.prefix_rank( 8, 'l') == 7);
        CHECK(vec.prefix_rank( 9, 'l') == 8);
        CHECK(vec.prefix_rank(10, 'l') == 8);

        CHECK(vec.prefix_rank( 0, 'o') == 0);
        CHECK(vec.prefix_rank( 1, 'o') == 1);
        CHECK(vec.prefix_rank( 2, 'o') == 2);
        CHECK(vec.prefix_rank( 3, 'o') == 3);
        CHECK(vec.prefix_rank( 4, 'o') == 4);
        CHECK(vec.prefix_rank( 5, 'o') == 5);
        CHECK(vec.prefix_rank( 6, 'o') == 6);
        CHECK(vec.prefix_rank( 7, 'o') == 7);
        CHECK(vec.prefix_rank( 8, 'o') == 8);
        CHECK(vec.prefix_rank( 9, 'o') == 9);
        CHECK(vec.prefix_rank(10, 'o') == 9);

        CHECK(vec.prefix_rank( 0, 't') ==  0);
        CHECK(vec.prefix_rank( 1, 't') ==  1);
        CHECK(vec.prefix_rank( 2, 't') ==  2);
        CHECK(vec.prefix_rank( 3, 't') ==  3);
        CHECK(vec.prefix_rank( 4, 't') ==  4);
        CHECK(vec.prefix_rank( 5, 't') ==  5);
        CHECK(vec.prefix_rank( 6, 't') ==  6);
        CHECK(vec.prefix_rank( 7, 't') ==  7);
        CHECK(vec.prefix_rank( 8, 't') ==  8);
        CHECK(vec.prefix_rank( 9, 't') ==  9);
        CHECK(vec.prefix_rank(10, 't') == 10);
    }

    SECTION("check all_ranks() is equal to prefix_rank() and rank()") {
        for (size_t idx{0}; idx < vec.size(); ++idx) {
            auto [rank, prefix] = vec.all_ranks_and_prefix_ranks(idx);
            for (size_t symb{1}; symb < Vector::Sigma; ++symb) {
                INFO(idx);
                INFO(symb);
                CHECK(rank[symb] == vec.rank(idx, symb));
                CHECK(prefix[symb] == vec.prefix_rank(idx, symb));
            }
        }
    }
}

TEMPLATE_TEST_CASE("check symbol vectors construction on text longer than 256 characters", "[SymbolVector]", ALLSYMBOLVECTORS) {
    using Vector = TestType;
    INFO(typeid(Vector).name());

    auto text = std::vector<uint8_t>{'H', 'a', 'l', 'l', 'o', ' ', 'W', 'e', 'l', 't',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    };

    auto vector = Vector{text};

    REQUIRE(vector.size() == text.size());

    SECTION("check that symbol() call works") {
        for (size_t i{0}; i < text.size(); ++i) {
            INFO(i);
            CHECK(vector.symbol(i) == text.at(i));
        }
    }
    CHECK(Vector::Sigma >= 128);

    SECTION("check all_ranks() is equal to prefix_rank() and rank()") {
        for (size_t idx{0}; idx < vector.size(); ++idx) {
            auto [rank, prefix] = vector.all_ranks_and_prefix_ranks(idx);
            for (size_t symb{1}; symb < Vector::Sigma; ++symb) {
                INFO(idx);
                INFO(symb);
                CHECK(rank[symb] == vector.rank(idx, symb));
                CHECK(prefix[symb] == vector.prefix_rank(idx, symb));
            }
        }
    }
}
