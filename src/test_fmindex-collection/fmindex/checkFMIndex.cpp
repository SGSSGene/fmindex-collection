// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "../string/utils.h"

#include <catch2/catch_all.hpp>
#include <fmindex-collection/fmindex/FMIndex.h>
#include <fmindex-collection/suffixarray/CSA.h>
#include <fmindex-collection/suffixarray/utils.h>

#include <fstream>

TEST_CASE("checking unidirectional fm index", "[fmindex]") {
    auto bwt    = std::vector<uint8_t>{'t', '\0', 'o', '\0', ' ', 'H', 'W', 'a', 'l', 'e', 'l', 'l'};
    auto sa     = std::vector<uint64_t>{ 10, 11, 5, 0,  6,  1,  7,  2,  3,  8,  4,  9 };

    SECTION("full sa") {
        auto bitStack = std::vector<bool>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            bitStack.push_back(true);
        }
        auto csa = fmc::CSA{sa, bitStack, /*.threadNbr=*/63, /*.seqCount=*/1};
        auto index = fmc::FMIndex<255>{bwt, fmc::suffixarray::convertCSAToAnnotatedDocument(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            auto [entry, offset] = index.locate(i);
            auto [seqId, seqPos] = entry;
            INFO(i);
            CHECK(seqId == 0);
            CHECK(seqPos == sa[i]-offset);
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

        auto csa = fmc::CSA{sa2, bitStack, /*.threadNbr=*/63, /*.seqCount=*/1};
        auto index = fmc::FMIndex<255>{bwt, fmc::suffixarray::convertCSAToAnnotatedDocument(csa)};

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

        auto csa = fmc::CSA{sa2, bitStack, /*.threadNbr=*/63, /*.seqCount=*/1};
        auto index = fmc::FMIndex<255>{bwt, fmc::suffixarray::convertCSAToAnnotatedDocument(csa)};

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

        auto csa = fmc::CSA{sa2, bitStack, /*.threadNbr=*/63, /*.seqCount=*/1};
        auto index = fmc::FMIndex<255>{bwt, fmc::suffixarray::convertCSAToAnnotatedDocument(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            auto [entry, offset] = index.locate(i);
            CHECK(entry == std::make_tuple(0, sa[i]-offset));
        }
    }

    SECTION("serialization/deserialization") {
        SECTION("serialize") {
            auto ofs = std::ofstream{"temp_test_serialization", std::ios::binary};
            auto bitStack = std::vector<bool>{};
            for (size_t i{0}; i < sa.size(); ++i) {
                bitStack.push_back(true);
            }
            auto csa = fmc::CSA{sa, bitStack, /*.threadNbr=*/63, /*.seqCount=*/1};
            auto index = fmc::FMIndex<255>{bwt, fmc::suffixarray::convertCSAToAnnotatedDocument(csa)};
            auto archive = cereal::BinaryOutputArchive{ofs};
            archive(index);
        }
        SECTION("deserialize") {
            auto ifs = std::ifstream{"temp_test_serialization", std::ios::binary};

            auto index = fmc::FMIndex<255>{};
            auto archive = cereal::BinaryInputArchive{ifs};
            archive(index);

            REQUIRE(index.size() == bwt.size());
            for (size_t i{0}; i < sa.size(); ++i) {
                auto [entry, offset] = index.locate(i);
                CHECK(entry == std::make_tuple(0, sa[i]-offset));
            }
        }
    }
}
