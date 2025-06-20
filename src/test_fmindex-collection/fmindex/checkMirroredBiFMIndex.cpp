// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "../string/utils.h"

#include <catch2/catch_all.hpp>
#include <fmindex-collection/fmindex/MirroredBiFMIndex.h>

TEST_CASE("checking mirrored bidirectional fm index", "[mirroredbifmindex]") {
    using String = fmindex_collection::string::PairedFlattenedBitvectors_512_64k<255>;

    // T = 011202110
    auto bwt    = std::vector<uint8_t>{1, 0, 2, 1, 2, 0, 1, 1, 0};
    auto sa     = std::vector<uint64_t>{ 7, 8, 3, 6, 5, 0, 1, 2, 4};

    SECTION("full sa") {
        auto bitStack = std::vector<bool>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            bitStack.push_back(true);
        }
        auto csa = fmindex_collection::CSA{sa, bitStack, /*.threadNbr=*/63, /*.seqCount=*/1};
        auto index = fmindex_collection::MirroredBiFMIndex<String>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            auto [entry, offset] = index.locate(i);
            CHECK(entry == std::make_tuple(0, sa[i]-offset));
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
        auto index = fmindex_collection::MirroredBiFMIndex<String>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            auto [entry, offset] = index.locate(i);
            CHECK(entry == std::make_tuple(0, sa[i]-offset));
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
            auto add = bool{i % 2 == 1};
            bitStack.push_back(add);
            if (add) {
                sa2.push_back(sa[i]);
            }
        }

        auto csa = fmindex_collection::CSA{sa2, bitStack, /*.threadNbr=*/63, /*.seqCount=*/1};
        auto index = fmindex_collection::MirroredBiFMIndex<String>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            auto [entry, offset] = index.locate(i);
            CHECK(entry == std::make_tuple(0, sa[i]-offset));
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
        auto index = fmindex_collection::MirroredBiFMIndex<String>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            auto [entry, offset] = index.locate(i);
            CHECK(entry == std::make_tuple(0, sa[i]-offset));
        }
    }
}
