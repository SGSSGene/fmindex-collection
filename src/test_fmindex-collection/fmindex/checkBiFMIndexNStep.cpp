// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "../string/allStrings.h"
#include "../string/utils.h"

#include <catch2/catch_all.hpp>
#include <fmindex-collection/fmindex/BiFMIndexNStep.h>
#include <fmindex-collection/suffixarray/CSA.h>
#include <fstream>

TEST_CASE("checking bidirectional fm index with nstep capabilities", "[bifmindex_nstep]") {
    auto bwt          = std::vector<uint8_t>{2, 2, 2, 3, 1, 3, 1, 1, 0, 1};
    auto bwtRev       = std::vector<uint8_t>{3, 3, 1, 2, 2, 1, 1, 0, 1, 2};
    auto bwt_nstep    = std::vector<uint8_t>{11,  9, 9, 12, 6, 13, 6, 7, 2, 5};
    auto bwtRev_nstep = std::vector<uint8_t>{13, 14, 7,  9, 9,  5, 6, 3, 6, 8};
    auto sa           = std::vector<uint64_t>{9, 5, 3, 1, 6, 8, 4, 2, 0, 7};

    SECTION("full sa") {
        auto bitStack = std::vector<bool>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            bitStack.push_back(true);
        }
        auto csa = fmc::CSA{sa, bitStack, /*.threadNbr=*/63, /*.seqCount=*/1};
        auto index = fmc::BiFMIndexNStep<4>{bwt, bwtRev, bwt_nstep, bwtRev_nstep, fmc::suffixarray::convertCSAToAnnotatedDocument(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i], 0));
        }
        for (size_t i{0}; i < bwt.size(); ++i) {
            INFO(i);
            CHECK(index.bwt.symbol(i) == bwt[i]);
            CHECK(index.bwtRev.symbol(i) == bwtRev[i]);
            CHECK(index.bwt_nstep.symbol(i) == bwt_nstep[i]);
            CHECK(index.bwtRev_nstep.symbol(i) == bwtRev_nstep[i]);
        }
    }

    SECTION("check sa creates the same transformation") {
        auto text = std::vector<std::vector<uint8_t>>{
            { 3, 1, 2, 1, 2, 1, 1, 3, 2}
        };
        auto index = fmc::BiFMIndexNStep<4>{text, /*.samplingRate=*/1, /*.threadNbr=*/1};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i], 0));
        }
        // Complecated, since they kind of have to be reversed.
        for (size_t i{0}; i < bwt.size(); ++i) {
            INFO(i);
            CHECK(index.bwt.symbol(i) == bwt[i]);
            CHECK(index.bwtRev.symbol(i) == bwtRev[i]);
            CHECK(index.bwt_nstep.symbol(i) == static_cast<size_t>(bwt_nstep[i]));
            CHECK(index.bwtRev_nstep.symbol(i) == static_cast<size_t>(bwtRev_nstep[i]));
        }
    }

#if 1

    SECTION("serialization/deserialization") {
        SECTION("serialize") {
            auto ofs = std::ofstream{"temp_test_serialization", std::ios::binary};

            auto bitStack = std::vector<bool>{};
            for (size_t i{0}; i < sa.size(); ++i) {
                bitStack.push_back(true);
            }
            auto csa = fmc::CSA{sa, bitStack, /*.threadNbr=*/63, /*.seqCount=*/1};
            auto index = fmc::BiFMIndexNStep<4>{bwt, bwtRev, bwt_nstep, bwtRev_nstep, fmc::suffixarray::convertCSAToAnnotatedDocument(csa)};
            auto archive = cereal::BinaryOutputArchive{ofs};
            archive(index);
        }
        SECTION("deserialize") {
            auto ifs = std::ifstream{"temp_test_serialization", std::ios::binary};

            auto index = fmc::BiFMIndexNStep<4>{};
            auto archive = cereal::BinaryInputArchive{ifs};
            archive(index);

            REQUIRE(index.size() == bwt.size());
            for (size_t i{0}; i < sa.size(); ++i) {
                CHECK(index.locate(i) == std::make_tuple(0, sa[i], 0));
            }
            for (size_t i{0}; i < bwt.size(); ++i) {
                CHECK(index.bwt.symbol(i) == bwt[i]);
                CHECK(index.bwtRev.symbol(i) == bwtRev[i]);
                CHECK(index.bwt_nstep.symbol(i) == bwt_nstep[i]);
                CHECK(index.bwtRev_nstep.symbol(i) == bwtRev_nstep[i]);
            }
        }
    }
    #endif
}
