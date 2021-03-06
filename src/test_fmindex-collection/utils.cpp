#include <fmindex-collection/utils.h>
#include <catch2/catch.hpp>

#include <iostream>


TEST_CASE("check creation of suffix array", "[createSA]") {
    auto input    = std::vector<uint8_t>{'H', 'a', 'l', 'l', 'o', ' ', 'W', 'e', 'l', 't', '\0', '\0'};
    auto expected = std::vector<int64_t>{ 11, 10, 5, 0,  6,  1,  7,  2,  3,  8,  4,  9 };

    auto output = fmindex_collection::createSA(input);
    CHECK(output == expected);
}

TEST_CASE("check creation of bwt", "[createBWT]") {
    auto text     = std::vector<uint8_t>{'H', 'a', 'l', 'l', 'o', ' ', 'W', 'e', 'l', 't', '\0', '\0'};
    auto expected = std::vector<uint8_t>{'\0', 't', 'o', '\0', ' ', 'H', 'W', 'a', 'l', 'e', 'l', 'l'};

    auto sa = fmindex_collection::createSA(text);
    auto output = fmindex_collection::createBWT(text, sa);

    CHECK(output == expected);
}
