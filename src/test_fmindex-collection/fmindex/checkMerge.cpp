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

    auto index1 = Index{data1, /*.samplingRate =*/ 1, /*.threadNbr =*/ 1};
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
        {0, 8},
        {1, 8},
        {0, 0},
        {0, 1},
        {0, 2},
        {0, 3},
        {0, 7},
        {1, 0},
        {1, 2},
        {1, 4},
        {1, 6},
        {0, 6},
        {0, 5},
        {0, 4},
        {1, 1},
        {1, 5},
        {1, 4},
        {1, 3},
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
    CHECK(texts[0] == std::vector<uint8_t>{0, 2, 1, 2, 1, 2, 1, 2, 2});
    CHECK(texts[1] == std::vector<uint8_t>{0, 1, 1, 1, 1, 2, 2, 2, 2});
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

    std::cout << "idx: " << "\n";
    size_t idx1{};
    for (auto _ : data1[0]) {
        (void)_;
        char c = index1.occ.symbol(idx1);
        idx1 = index1.occ.rank(idx1, c);
        std::cout << (int)c << " " << idx1 << "\n";
    }

    for (size_t i{0}; i < index1.size(); ++i) {
        auto idx0 = index1.occ.rank(i, 0);
        auto idx1 = index1.occ.rank(i, 1);
        auto idx2 = index1.occ.rank(i, 2);

        auto [seq, pos] = index1.locate(i);

        std::cout << idx0 << " " << idx1 << " " << idx2 << " : " << seq << " " << pos << "\n";
    }

    auto index12 = merge(index1, index2);
    std::cout << "reconstruct 12:\n";
    for (auto const& t : reconstructText(index12)) {
        for (auto c : t) {
            std::cout << (int)c;
        }
        std::cout << "\n";
    }
    for (size_t i{0}; i < index12.size(); ++i) {
        auto idx0 = index12.occ.rank(i, 0);
        auto idx1 = index12.occ.rank(i, 1);
        auto idx2 = index12.occ.rank(i, 2);

        auto [seq, pos] = index12.locate(i);

        std::cout << idx0 << " " << idx1 << " " << idx2 << " : " << seq << " " << pos << "\n";
    }



    auto index123 = merge(index12, index3);


    std::cout << "reconstruct 123:\n";
    for (auto const& t : reconstructText(index123)) {
        for (auto c : t) {
            std::cout << (int)c;
        }
        std::cout << "\n";
    }
    for (size_t i{0}; i < index123.size(); ++i) {
        auto idx0 = index123.occ.rank(i, 0);
        auto idx1 = index123.occ.rank(i, 1);
        auto idx2 = index123.occ.rank(i, 2);

        auto [seq, pos] = index123.locate(i);

        std::cout << idx0 << " " << idx1 << " " << idx2 << " : " << seq << " " << pos << "\n";
    }


}
