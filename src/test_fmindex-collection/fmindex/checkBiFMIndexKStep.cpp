// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "../string/allStrings.h"
#include "../string/utils.h"

#include <catch2/catch_all.hpp>
#include <fmindex-collection/fmindex/BiFMIndexKStep.h>
#include <fmindex-collection/suffixarray/CSA.h>
#include <fstream>

    // $ 31212113 2    $3 121211 32   -> 11
    // 1 132$3121 2    11 32$312 12       9
    // 1 21132$31 2    12 1132$3 12       9
    // 1 2121132$ 3    12 121132 $3      12
    // 1 32$31212 1    13 2$3121 21       6
    // 2 $3121211 3    2$ 312121 13      13
    // 2 1132$312 1    21 132$31 21       6
    // 2 121132$3 1    21 21132$ 31       7
    // 3 12121132 $    31 212113 2$       2
    // 3 2$312121 1    32 $31212 11       5
    //
    //
    // $ 23112121 3    $2 311212 13  -> 13
    // 1 121213$2 3    11 21213$ 23     14
    // 1 21213$23 1    12 1213$2 31      7
    // 1 213$2311 2    12 13$231 12      9
    // 1 3$231121 2    13 $23112 12      9
    // 2 1213$231 1    21 213$23 11      5
    // 2 13$23112 1    21 3$2311 21      6
    // 2 31121213 $    23 112121 3$      3
    // 3 $2311212 1    3$ 231121 21      6
    // 3 1121213$ 2    31 121213 $2      8



TEST_CASE("checking bidirectional fm index with kstep capabilities - ctor", "[bifmindex_kstep]") {
    auto bwt              = std::vector<uint8_t>{2, 2, 2, 3, 1, 3, 1, 1, 0, 1};
    auto bwtRev           = std::vector<uint8_t>{3, 3, 1, 2, 2, 1, 1, 0, 1, 2};
    auto bwt_kstep        = std::vector<uint8_t>{11,  9, 9, 12, 6, 13, 6, 7, 2, 5};
    auto bwtRev_kstep     = std::vector<uint8_t>{13, 14, 7,  9, 9,  5, 6, 3, 6, 8};
    auto sa               = std::vector<uint64_t>{9, 5, 3, 1, 6, 8, 4, 2, 0, 7};
    auto input            = std::vector<std::vector<uint8_t>>{std::vector<uint8_t>{3, 1, 2, 1, 2, 1, 1, 3, 2}};
//    auto at               = std::vector<uint64_t>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
//    auto isKStep          = std::vector<bool>{true, false, true, false, true, false, true, false, true, false};

//    BiFMIndexKStep(Sequence auto const& _sequence, SparseArray const& _annotatedSequence, auto _annotatedSequenceIsKStep, size_t _threadNbr, bool includeReversedInput = false) {
//    auto index = fmc::BiFMIndexKStep<4>{input, at, isKStep, 1, false};
//    auto index = fmc::BiFMIndexKStep<4>(Sequences auto const& _input, size_t samplingRate, size_t threadNbr, size_t seqOffset = 0, bool includeReversedInput = false) {
//    fmc::string::AdapterStringKStep<16, 2, fmc::string::InterleavedBitvectorPrefix16>{input2};/
    using Index = fmc::BiFMIndexKStep<4, fmc::string::InterleavedBitvectorPrefix16>;
//    using Index = fmc::BiFMIndexKStep<4>;
    auto index = Index{input, 16, 1};
    for (size_t i{0}; i < bwt.size(); ++i) {
        INFO(i);
        CHECK(index.bwt_kstep.symbol_limit<2>(i) == bwt[i]);
        CHECK(index.bwtRev_kstep.symbol_limit<2>(i) == bwtRev[i]);
        CHECK(index.bwt_kstep.symbol(i) == bwt_kstep[i]);
        CHECK(index.bwtRev_kstep.symbol(i) == bwtRev_kstep[i]);
        auto [seqId, pos, offset] = index.locate(i);
        CHECK(seqId == 0);
        CHECK(sa[i] == pos + offset);

    }

    mmser::saveFile("/tmp/test.mmser", index);
    auto [index2, storageManager] = mmser::loadFile<Index>("/tmp/test.mmser");
    for (size_t i{0}; i < bwt.size(); ++i) {
        INFO(i);
        CHECK(index2.bwt_kstep.symbol_limit<2>(i) == bwt[i]);
        CHECK(index2.bwtRev_kstep.symbol_limit<2>(i) == bwtRev[i]);
        CHECK(index2.bwt_kstep.symbol(i) == bwt_kstep[i]);
        CHECK(index2.bwtRev_kstep.symbol(i) == bwtRev_kstep[i]);
        auto [seqId, pos, offset] = index2.locate(i);
        CHECK(seqId == 0);
        CHECK(sa[i] == pos + offset);
    }
}

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
            CHECK(index.bwt_kstep.symbol_limit<2>(i) == bwt[i]);
            CHECK(index.bwtRev_kstep.symbol_limit<2>(i) == bwtRev[i]);
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
            CHECK(index.bwt_kstep.symbol_limit<2>(i) == bwt[i]);
            CHECK(index.bwtRev_kstep.symbol_limit<2>(i) == bwtRev[i]);
            CHECK(index.bwt_kstep.symbol(i) == static_cast<size_t>(bwt_kstep[i]));
            CHECK(index.bwtRev_kstep.symbol(i) == static_cast<size_t>(bwtRev_kstep[i]));
        }
    }
    SECTION("check creation for InterleavedBitvectorPrefix16") {
        auto text = std::vector<std::vector<uint8_t>>{
            { 3, 1, 2, 1, 2, 1, 1, 3, 2}
        };
        auto index = fmc::BiFMIndexKStep<4, fmc::string::InterleavedBitvectorPrefix16>{text, /*.samplingRate=*/1, /*.threadNbr=*/1};
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
