// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <catch2/catch_all.hpp>
#include <fmindex-collection/suffixarray/CSA.h>
#include <fmindex-collection/suffixarray/DenseCSA.h>
#include <fmindex-collection/utils.h>

TEST_CASE("check Compressed Suffix Array (CSA)", "[suffixarray][CSA]") {
    // SA of 'Hello$World$':
    //               S P  SA
    // $Hello$World  1 5  11
    // $World$Hello  0 5   5
    // Hello$World$  0 0   0
    // World$Hello$  1 0   6
    // d$Hello$Worl  1 4  10
    // ello$World$H  0 1   1
    // ld$Hello$Wor  1 3   9
    // llo$World$He  0 2   2
    // lo$World$Hel  0 3   3
    // o$World$Hell  0 4   4
    // orld$Hello$W  1 1   7
    // rld$Hello$Wo  1 2   8


    auto [totalSize, input, inputSizes] = fmindex_collection::createSequences(fmindex_collection::test::createInputData({"Hello", "World"}));
    auto sa                             = fmindex_collection::createSA64(input, 1);


    auto check = [&](auto const& csa, auto const& expected) {
        REQUIRE(sa.size() == expected.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            REQUIRE(csa.value(i) == expected[i]);
            if (expected[i]) {
                CHECK(*csa.value(i) == *expected[i]);
            }
        }
    };

    SECTION("sampling 3") {
        auto csa = fmindex_collection::CSA {sa, /*.samplingRate =*/ 3, inputSizes};

        auto expected = std::vector<std::optional<std::tuple<size_t, size_t>>>(12);
        expected[2] = {0, 0};
        expected[3] = {1, 0};
        expected[6] = {1, 3};
        expected[8] = {0, 3};
        check(csa, expected);
    }

    SECTION("sampling 4") {
        auto csa = fmindex_collection::CSA {sa, /*.samplingRate =*/ 4, inputSizes};

        auto expected = std::vector<std::optional<std::tuple<size_t, size_t>>>(12);
        expected[2] = {0, 0};
        expected[3] = {1, 0};
        expected[4] = {1, 4};
        expected[9] = {0, 4};
        check(csa, expected);
    }

    SECTION("sampling 5") {
        auto csa = fmindex_collection::CSA {sa, /*.samplingRate =*/ 5, inputSizes};

        auto expected = std::vector<std::optional<std::tuple<size_t, size_t>>>(12);
        expected[0] = {1, 5};
        expected[1] = {0, 5};
        expected[2] = {0, 0};
        expected[3] = {1, 0};
        check(csa, expected);
    }

    SECTION("sampling 8") {
        auto csa = fmindex_collection::CSA {sa, /*.samplingRate =*/ 8, inputSizes};

        auto expected = std::vector<std::optional<std::tuple<size_t, size_t>>>(12);
        expected[2] = {0, 0};
        expected[3] = {1, 0};
        check(csa, expected);
    }
}

TEST_CASE("check Dense Compressed Suffix Array (DenseCSA)", "[suffixarray][DenseCSA]") {
    // SA of 'Hello$World$':
    //               S P  SA
    // $Hello$World  1 5  11
    // $World$Hello  0 5   5
    // Hello$World$  0 0   0
    // World$Hello$  1 0   6
    // d$Hello$Worl  1 4  10
    // ello$World$H  0 1   1
    // ld$Hello$Wor  1 3   9
    // llo$World$He  0 2   2
    // lo$World$Hel  0 3   3
    // o$World$Hell  0 4   4
    // orld$Hello$W  1 1   7
    // rld$Hello$Wo  1 2   8


    auto [totalSize, input, inputSizes] = fmindex_collection::createSequences(fmindex_collection::test::createInputData({"Hello", "World"}));
    auto sa                             = fmindex_collection::createSA64(input, 1);


    auto check = [&](auto const& csa, auto const& expected) {
        REQUIRE(sa.size() == expected.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            REQUIRE(csa.value(i) == expected[i]);
            if (expected[i]) {
                CHECK(*csa.value(i) == *expected[i]);
            }
        }
    };

    SECTION("sampling 3") {
        auto csa = fmindex_collection::DenseCSA {sa, /*.samplingRate =*/ 3, inputSizes};

        auto expected = std::vector<std::optional<std::tuple<size_t, size_t>>>(12);
        expected[2] = {0, 0};
        expected[3] = {1, 0};
        expected[6] = {1, 3};
        expected[8] = {0, 3};
        check(csa, expected);
    }

    SECTION("sampling 4") {
        auto csa = fmindex_collection::DenseCSA {sa, /*.samplingRate =*/ 4, inputSizes};

        auto expected = std::vector<std::optional<std::tuple<size_t, size_t>>>(12);
        expected[2] = {0, 0};
        expected[3] = {1, 0};
        expected[4] = {1, 4};
        expected[9] = {0, 4};
        check(csa, expected);
    }

    SECTION("sampling 5") {
        auto csa = fmindex_collection::DenseCSA {sa, /*.samplingRate =*/ 5, inputSizes};

        auto expected = std::vector<std::optional<std::tuple<size_t, size_t>>>(12);
        expected[0] = {1, 5};
        expected[1] = {0, 5};
        expected[2] = {0, 0};
        expected[3] = {1, 0};
        check(csa, expected);
    }

    SECTION("sampling 8") {
        auto csa = fmindex_collection::DenseCSA {sa, /*.samplingRate =*/ 8, inputSizes};

        auto expected = std::vector<std::optional<std::tuple<size_t, size_t>>>(12);
        expected[2] = {0, 0};
        expected[3] = {1, 0};
        check(csa, expected);
    }

}
