// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include <catch2/catch_all.hpp>
#include <fmindex-collection/utils.h>

#include <iostream>

TEST_CASE("check creation of suffix array", "[createSA]") {
    auto input    = std::vector<uint8_t>{'H', 'a', 'l', 'l', 'o', ' ', 'W', 'e', 'l', 't', '\0', '\0'};
    auto expected = std::vector<uint64_t>{ 11, 10, 5, 0,  6,  1,  7,  2,  3,  8,  4,  9 };

    auto output = fmindex_collection::createSA(input, 1);
    CHECK(output == expected);
}

TEST_CASE("check creation of bwt", "[createBWT]") {
    auto text     = std::vector<uint8_t>{'H', 'a', 'l', 'l', 'o', ' ', 'W', 'e', 'l', 't', '\0', '\0'};
    auto expected = std::vector<uint8_t>{'\0', 't', 'o', '\0', ' ', 'H', 'W', 'a', 'l', 'e', 'l', 'l'};

    auto sa = fmindex_collection::createSA(text, 1);
    auto output = fmindex_collection::createBWT(text, sa);

    CHECK(output == expected);
}
