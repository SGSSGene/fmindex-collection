// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include "../bitvector/allBitVectors.h"

#include <catch2/catch_all.hpp>
#include <fmindex-collection/fmindex/BinaryMirroredBiFMIndex.h>

TEMPLATE_TEST_CASE("checking binary mirrored bidirectional fm index", "[BinaryMirroredBiFMIndex]", ALLBITVECTORS) {
    using String = TestType;
    // T = 011101010
    auto bwt    = std::vector<uint8_t> {1,  1,  1,  1,  0,  1,  1,  0,  0,  1,  0,  1,  0,  0,  1,  1,  0,  0};
    auto sa     = std::vector<uint64_t>{8, 17,  6,  4,  9, 11, 13,  0,  7, 16,  5,  3, 10, 12, 15,  2, 14,  1};
    assert(bwt.size() == sa.size());

    SECTION("full sa") {
        auto bitStack = std::vector<bool>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            bitStack.push_back(true);
        }
        auto csa = fmindex_collection::CSA{sa, bitStack, /*.threadNbr=*/63, /*.seqCount=*/1};
        auto index = fmindex_collection::BinaryMirroredBiFMIndex<String>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }

    SECTION("sa with only every second value given - sa sampled") {
        auto bitStack = std::vector<bool>{};
        auto sa2 = std::vector<uint64_t>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            auto add = bool{i % 2 == 0} || (sa[i] == 0);
            bitStack.push_back(add);
            if (add) {
                sa2.push_back(sa[i]);
            }
        }

        auto csa = fmindex_collection::CSA{sa2, bitStack, /*.threadNbr=*/63, /*.seqCount=*/1};
        auto index = fmindex_collection::BinaryMirroredBiFMIndex<String>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            if (index.locate(i) != std::make_tuple(0, sa[i])) {
                index.locate(i);
            }
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
            auto res = index.single_locate_step(i);
            INFO(i);
            INFO(sa[i]);
            if (sa[i] == 0 || i % 2 == 0) {
                REQUIRE(res);
                CHECK(*res == std::make_tuple(0, sa[i]));
            } else {
                CHECK(!res);
            }
        }
    }

    SECTION("sa with only every second value given - sa sampled - uneven") {
        auto bitStack = std::vector<bool>{};
        auto sa2 = std::vector<uint64_t>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            auto add = bool{i % 2 == 1} || (sa[i] == 0);
            bitStack.push_back(add);
            if (add) {
                sa2.push_back(sa[i]);
            }
        }

        auto csa = fmindex_collection::CSA{sa2, bitStack, /*.threadNbr=*/63, /*.seqCount=*/1};
        auto index = fmindex_collection::BinaryMirroredBiFMIndex<String>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            if (index.locate(i) != std::make_tuple(0, sa[i])) {
                index.locate(i);
            }
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }


    SECTION("sa with only every second value given - text sampled") {
        auto bitStack = std::vector<bool>{};
        auto sa2 = std::vector<uint64_t>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            auto add = bool{sa[i] % 2 == 0};
            bitStack.push_back(add);
            if (add) {
                sa2.push_back(sa[i]);
            }
        }

        auto csa = fmindex_collection::CSA{sa2, bitStack, /*.threadNbr=*/63, /*.seqCount=*/1};
        auto index = fmindex_collection::BinaryMirroredBiFMIndex<String>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }

    }
}
