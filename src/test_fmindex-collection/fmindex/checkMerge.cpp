// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include "../occtables/allTables.h"

#include <catch2/catch_all.hpp>
#include <fmindex-collection/fmindex/FMIndex.h>
#include <fmindex-collection/fmindex/merge.h>
#include <fmindex-collection/suffixarray/DenseCSA.h>

TEST_CASE("checking merging of fmindices", "[FMIndex][merge]") {

    auto data1 = std::vector<std::vector<uint8_t>>{std::vector<uint8_t>{1, 1, 1, 1, 2, 2, 2, 2}};
    auto data2 = std::vector<std::vector<uint8_t>>{std::vector<uint8_t>{2, 1, 2, 1, 2, 1, 2, 2}};

    using OccTable = fmindex_collection::occtable::Bitvector<3>;
    using Index = fmindex_collection::FMIndex<OccTable>;

    auto index1 = Index{data1, /*.samplingRate =*/ 2, /*.threadNbr =*/ 1};
    auto index2 = Index{data2, /*.samplingRate =*/ 2, /*.threadNbr =*/ 1};

    auto index12 = merge(index1, index2);

    auto expectedRanks = std::vector<std::tuple<size_t, size_t, size_t>> {
        { 0, 2,  9},
        { 0, 2, 10},
        { 0, 2, 11},
        { 1, 2, 11},
        { 1, 3, 11},
        { 1, 4, 11},
        { 1, 4, 12},
        { 1, 4, 13},
        { 1, 4, 14},
        { 1, 5, 14},
        { 1, 5, 15},
        { 1, 5, 16},
        { 2, 5, 16},
        { 2, 6, 16},
        { 2, 7, 16},
        { 2, 8, 16},
        { 2, 8, 17},
        { 2, 8, 18},
    };
    auto expectedSA = std::vector<std::tuple<size_t, size_t>> {
        {1, 8},
        {0, 8},
        {0, 0},
        {0, 1},
        {0, 2},
        {1, 1},
        {1, 3},
        {1, 5},
        {0, 3},
        {1, 7},
        {0, 7},
        {1, 0},
        {1, 2},
        {1, 4},
        {1, 6},
        {0, 6},
        {0, 5},
        {0, 4},
    };
    for (size_t i{0}; i < index12.size(); ++i) {
        auto idx0 = index12.occ.rank(i, 0);
        auto idx1 = index12.occ.rank(i, 1);
        auto idx2 = index12.occ.rank(i, 2);
        auto [seq, pos] = index12.locate(i);
        INFO(i);
        CHECK(std::get<0>(expectedRanks[i]) == idx0);
        CHECK(std::get<1>(expectedRanks[i]) == idx1);
        CHECK(std::get<2>(expectedRanks[i]) == idx2);
        CHECK(std::get<0>(expectedSA[i]) == seq);
        CHECK(std::get<1>(expectedSA[i]) == pos);
    }


    auto texts = reconstructText(index12);
    REQUIRE(texts.size() == 2);
    CHECK(texts[0] == data1[0]);
    CHECK(texts[1] == data2[0]);
}

TEST_CASE("checking merging of fmindices", "[BiFMIndex][merge]") {

    auto data1 = std::vector<std::vector<uint8_t>>{std::vector<uint8_t>{1, 1, 1, 1, 2, 2, 2, 2}};
    auto data2 = std::vector<std::vector<uint8_t>>{std::vector<uint8_t>{2, 1, 2, 1, 2, 1, 2, 2}};
    auto data3 = std::vector<std::vector<uint8_t>>{std::vector<uint8_t>{1, 1, 2, 2, 1, 1, 2, 2}};

    using OccTable = fmindex_collection::occtable::Bitvector<3>;
    using Index = fmindex_collection::BiFMIndex<OccTable, fmindex_collection::DenseCSA>;

    auto index1 = Index{data1, /*.samplingRate =*/ 2, /*.threadNbr =*/ 1};
    auto index2 = Index{data2, /*.samplingRate =*/ 2, /*.threadNbr =*/ 1};
    auto index3 = Index{data3, /*.samplingRate =*/ 2, /*.threadNbr =*/ 1};

    SECTION("merging index1 and index2 into index12") {
        auto index12 = merge(index1, index2);

        auto expectedRanks = std::vector<std::tuple<size_t, size_t, size_t>> {
            { 0, 2,  9},
            { 0, 2, 10},
            { 0, 2, 11},
            { 1, 2, 11},
            { 1, 3, 11},
            { 1, 4, 11},
            { 1, 4, 12},
            { 1, 4, 13},
            { 1, 4, 14},
            { 1, 5, 14},
            { 1, 5, 15},
            { 1, 5, 16},
            { 2, 5, 16},
            { 2, 6, 16},
            { 2, 7, 16},
            { 2, 8, 16},
            { 2, 8, 17},
            { 2, 8, 18},
        };
        auto expectedSA = std::vector<std::tuple<size_t, size_t>> {
            {1, 8},
            {0, 8},
            {0, 0},
            {0, 1},
            {0, 2},
            {1, 1},
            {1, 3},
            {1, 5},
            {0, 3},
            {1, 7},
            {0, 7},
            {1, 0},
            {1, 2},
            {1, 4},
            {1, 6},
            {0, 6},
            {0, 5},
            {0, 4},
        };
        for (size_t i{0}; i < index12.size(); ++i) {
            auto idx0 = index12.occ.rank(i, 0);
            auto idx1 = index12.occ.rank(i, 1);
            auto idx2 = index12.occ.rank(i, 2);
            auto [seq, pos] = index12.locate(i);
            INFO(i);
            CHECK(std::get<0>(expectedRanks[i]) == idx0);
            CHECK(std::get<1>(expectedRanks[i]) == idx1);
            CHECK(std::get<2>(expectedRanks[i]) == idx2);
            CHECK(std::get<0>(expectedSA[i]) == seq);
            CHECK(std::get<1>(expectedSA[i]) == pos);
        }


        auto texts = reconstructText(index12);
        REQUIRE(texts.size() == 2);
        CHECK(texts[0] == data1[0]);
        CHECK(texts[1] == data2[0]);

        SECTION("merging index12 and index3 into index123") {
            auto index123 = merge(index12, index3);

            auto expectedRanks = std::vector<std::tuple<size_t, size_t, size_t>> {
                {0,  3, 14},
                {0,  3, 15},
                {0,  3, 16},
                {0,  3, 17},
                {1,  3, 17},
                {1,  4, 17},
                {1,  4, 18},
                {2,  4, 18},
                {2,  5, 18},
                {2,  5, 19},
                {2,  5, 20},
                {2,  6, 20},
                {2,  6, 21},
                {2,  7, 21},
                {2,  8, 21},
                {2,  8, 22},
                {2,  8, 23},
                {2,  8, 24},
                {2,  8, 25},
                {3,  8, 25},
                {3,  9, 25},
                {3, 10, 25},
                {3, 11, 25},
                {3, 12, 25},
                {3, 12, 26},
                {3, 13, 26},
                {3, 13, 27},
            };
            auto expectedSA = std::vector<std::tuple<size_t, size_t>> {
                {2, 8},
                {1, 8},
                {0, 8},
                {0, 0},
                {0, 1},
                {2, 4},
                {2, 0},
                {0, 2},
                {1, 1},
                {1, 3},
                {2, 5},
                {1, 5},
                {2, 1},
                {0, 3},
                {2, 7},
                {1, 7},
                {0, 7},
                {2, 3},
                {1, 0},
                {1, 2},
                {1, 4},
                {2, 6},
                {1, 6},
                {0, 6},
                {2, 2},
                {0, 5},
                {0, 4},
            };
            for (size_t i{0}; i < index123.size(); ++i) {
                auto idx0 = index123.occ.rank(i, 0);
                auto idx1 = index123.occ.rank(i, 1);
                auto idx2 = index123.occ.rank(i, 2);
                auto [seq, pos] = index123.locate(i);
                INFO(i);
                CHECK(std::get<0>(expectedRanks[i]) == idx0);
                CHECK(std::get<1>(expectedRanks[i]) == idx1);
                CHECK(std::get<2>(expectedRanks[i]) == idx2);
                CHECK(std::get<0>(expectedSA[i]) == seq);
                CHECK(std::get<1>(expectedSA[i]) == pos);
            }


            auto texts = reconstructText(index123);
            REQUIRE(texts.size() == 3);
            CHECK(texts[0] == data1[0]);
            CHECK(texts[1] == data2[0]);
            CHECK(texts[2] == data3[0]);
        }
    }
}
