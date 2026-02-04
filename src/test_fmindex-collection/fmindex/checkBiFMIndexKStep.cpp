// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "../string/allStrings.h"
#include "../string/utils.h"

#include <catch2/catch_all.hpp>
#include <fmindex-collection/fmindex/BiFMIndexKStep.h>
#include <fmindex-collection/suffixarray/CSA.h>
#include <fstream>

TEST_CASE("checking bidirectional fm index with kstep capabilities", "[bifmindex_kstep]") {
    auto bwt              = std::vector<uint8_t>{2, 2, 2, 3, 1, 3, 1, 1, 0, 1};
    auto bwtRev           = std::vector<uint8_t>{3, 3, 1, 2, 2, 1, 1, 0, 1, 2};
    auto bwt_kstep        = std::vector<uint8_t>{11,  9, 9, 12, 6, 13, 6, 7, 2, 5};
    auto bwtRev_kstep     = std::vector<uint8_t>{13, 14, 7,  9, 9,  5, 6, 3, 6, 8};
    auto sa               = std::vector<uint64_t>{9, 5, 3, 1, 6, 8, 4, 2, 0, 7};
    auto isKStepAnnotated = fmc::VectorBool{};
    isKStepAnnotated.resize(10, true);

    SECTION("full sa") {
        auto bitStack = std::vector<bool>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            bitStack.push_back(true);
        }
        auto csa = fmc::CSA{sa, bitStack, /*.threadNbr=*/63, /*.seqCount=*/1};
        auto index = fmc::BiFMIndexKStep<4>{bwt, bwt_kstep, bwtRev_kstep, fmc::suffixarray::convertCSAToAnnotatedDocument(csa), isKStepAnnotated};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i], 0));
        }
        for (size_t i{0}; i < bwt.size(); ++i) {
            INFO(i);
            //CHECK(index.bwt.symbol(i) == bwt[i]);
            //CHECK(index.bwtRev.symbol(i) == bwtRev[i]);
            CHECK(index.bwt_kstep.symbol(i) == bwt_kstep[i]);
            CHECK(index.bwtRev_kstep.symbol(i) == bwtRev_kstep[i]);
        }
    }

    SECTION("check sa creates the same transformation") {
        auto text = std::vector<std::vector<uint8_t>>{
            { 3, 1, 2, 1, 2, 1, 1, 3, 2}
        };
        auto index = fmc::BiFMIndexKStep<4>{text, /*.samplingRate=*/1, /*.threadNbr=*/1};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i], 0));
        }
        // Complicated, since they kind of have to be reversed.
        for (size_t i{0}; i < bwt.size(); ++i) {
            INFO(i);
//            CHECK(index.bwt.symbol(i) == bwt[i]);
//            CHECK(index.bwtRev.symbol(i) == bwtRev[i]);
            CHECK(index.bwt_kstep.symbol(i) == static_cast<size_t>(bwt_kstep[i]));
            CHECK(index.bwtRev_kstep.symbol(i) == static_cast<size_t>(bwtRev_kstep[i]));
        }
    }

    SECTION("sa with only every third value given - text sampled") {
        auto text = std::vector<std::vector<uint8_t>>{
            { 3, 1, 2, 1, 2, 1, 1, 3, 2}
        };
        auto index = fmc::BiFMIndexKStep<4>{text, /*.samplingRate=*/5, /*.threadNbr=*/1};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            auto [seqId, pos, offset] = index.locate(i);
            INFO("i " << i);
            INFO("pos " << pos);
            INFO("offset " << offset);
            CHECK(seqId == 0);
            CHECK(pos+offset == sa[i]);
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
            auto index = fmc::BiFMIndexKStep<4>{bwt, bwt_kstep, bwtRev_kstep, fmc::suffixarray::convertCSAToAnnotatedDocument(csa), isKStepAnnotated};
            auto archive = cereal::BinaryOutputArchive{ofs};
            archive(index);
        }
        SECTION("deserialize") {
            auto ifs = std::ifstream{"temp_test_serialization", std::ios::binary};

            auto index = fmc::BiFMIndexKStep<4>{};
            auto archive = cereal::BinaryInputArchive{ifs};
            archive(index);

            REQUIRE(index.size() == bwt.size());
            for (size_t i{0}; i < sa.size(); ++i) {
                CHECK(index.locate(i) == std::make_tuple(0, sa[i], 0));
            }
            for (size_t i{0}; i < bwt.size(); ++i) {
                CHECK(index.bwt_kstep.symbol(i) == bwt_kstep[i]);
                CHECK(index.bwtRev_kstep.symbol(i) == bwtRev_kstep[i]);
            }
        }
    }
}
