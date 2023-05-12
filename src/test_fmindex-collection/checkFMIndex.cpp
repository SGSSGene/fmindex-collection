// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#include "allTables.h"

#include <fmindex-collection/FMIndex.h>
#include <catch2/catch.hpp>

TEMPLATE_TEST_CASE("checking unidirectional fm index", "[FMIndex]", ALLTABLES) {
    using OccTable = TestType;

    auto bwt    = std::vector<uint8_t>{'t', '\0', 'o', '\0', ' ', 'H', 'W', 'a', 'l', 'e', 'l', 'l'};
    auto sa     = std::vector<uint64_t>{ 10, 11, 5, 0,  6,  1,  7,  2,  3,  8,  4,  9 };

    SECTION("full sa") {
        auto bitStack = fmindex_collection::BitStack{};
        for (size_t i{0}; i < sa.size(); ++i) {
            bitStack.push(true);
        }
        auto csa = fmindex_collection::CSA{sa, bitStack, 1, 63};
        auto index = fmindex_collection::FMIndex<OccTable>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            auto [seqId, seqPos] = index.locate(i);
            INFO(i);
            CHECK(seqId == 0);
            CHECK(seqPos == sa[i]);
        }
    }

    SECTION("sa with only every second value given - sa sampled") {
        auto bitStack = fmindex_collection::BitStack{};
        auto sa2 = std::vector<uint64_t>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            auto add = bool{i % 2 == 0} || (sa[i] == 0);
            bitStack.push(add);
            if (add) {
                sa2.push_back(sa[i]);
            }
        }

        auto csa = fmindex_collection::CSA{sa2, bitStack, 2, 63};
        auto index = fmindex_collection::FMIndex<OccTable>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }

    SECTION("sa with only every second value given - sa sampled - uneven") {
        auto bitStack = fmindex_collection::BitStack{};
        auto sa2 = std::vector<uint64_t>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            auto add = bool{i % 2 == 1};
            bitStack.push(add);
            if (add) {
                sa2.push_back(sa[i]);
            }
        }

        auto csa = fmindex_collection::CSA{sa2, bitStack, 2, 63};
        auto index = fmindex_collection::FMIndex<OccTable>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }


    SECTION("sa with only every second value given - text sampled") {
        auto bitStack = fmindex_collection::BitStack{};
        auto sa2 = std::vector<uint64_t>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            auto add = bool{sa[i] % 2 == 0};
            bitStack.push(add);
            if (add) {
                sa2.push_back(sa[i]);
            }
        }

        auto csa = fmindex_collection::CSA{sa2, bitStack, 2, 63};
        auto index = fmindex_collection::FMIndex<OccTable>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }

    }
}
