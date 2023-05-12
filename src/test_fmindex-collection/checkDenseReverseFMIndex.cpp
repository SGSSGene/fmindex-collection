// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#include "allTables.h"

#include <fmindex-collection/ReverseFMIndex.h>
#include <fmindex-collection/DenseCSA.h>
#include <catch2/catch.hpp>

TEMPLATE_TEST_CASE("checking dense reverse fm index", "[DenseReverseFMIndex]", ALLTABLES) {
    using OccTable = TestType;
    using DenseVector = fmindex_collection::DenseVector;

    INFO("OccTable " << typeid(OccTable).name());


    auto bwt = std::vector<uint8_t>{'H', '\0', 'W', 'a', 'e', 'l', 'l', 'l', 't', 'o', ' ', '\0'};
    auto sa  = DenseVector{0, 11, 6, 1, 7, 2, 8, 3, 9, 4, 5, 10};

    SECTION("full sa") {
        auto bitStack = fmindex_collection::BitStack{};
        for (size_t i{0}; i < sa.size(); ++i) {
            bitStack.push(true);
        }
        auto csa = fmindex_collection::DenseCSA{sa, bitStack, /*samplingRate*/ 1};
        auto index = fmindex_collection::ReverseFMIndex<OccTable, fmindex_collection::DenseCSA>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }

    SECTION("sa with only every second value given - sa sampled") {
        auto bitStack = fmindex_collection::BitStack{};
        auto sa2 = DenseVector(sa.entry_size());
        for (size_t i{0}; i < sa.size(); ++i) {
            auto add = bool{i % 2 == 0} || (sa[i] == sa.size()-1);
            bitStack.push(add);
            if (add) {
                sa2.push_back(sa[i]);
            }
        }

        auto csa = fmindex_collection::DenseCSA{sa2, bitStack, /*samplingRate*/ 2};
        auto index = fmindex_collection::ReverseFMIndex<OccTable, fmindex_collection::DenseCSA>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            INFO(i);
            INFO(std::get<1>(index.locate(i)));
            INFO(sa[i]);
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }

    SECTION("sa with only every second value given - sa sampled - uneven") {
        auto bitStack = fmindex_collection::BitStack{};
        auto sa2 = DenseVector(sa.entry_size());
        for (size_t i{0}; i < sa.size(); ++i) {
            auto add = bool{i % 2 == 1} || (sa[i] == sa.size()-1);
            bitStack.push(add);
            if (add) {
                sa2.push_back(sa[i]);
            }
        }

        auto csa = fmindex_collection::DenseCSA{sa2, bitStack, /*samplingRate*/ 2};
        auto index = fmindex_collection::ReverseFMIndex<OccTable, fmindex_collection::DenseCSA>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            INFO(i);
            INFO(std::get<1>(index.locate(i)));
            INFO(sa[i]);
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }


    SECTION("sa with only every second value given - text sampled") {
        auto bitStack = fmindex_collection::BitStack{};
        auto sa2 = DenseVector(sa.entry_size());
        for (size_t i{0}; i < sa.size(); ++i) {
            auto add = bool{sa[i] % 2 == 0} || (sa[i] == 11);
            bitStack.push(add);
            if (add) {
                sa2.push_back(sa[i]);
            }
        }

        auto csa = fmindex_collection::DenseCSA{sa2, bitStack, /*samplingRate*/ 2};
        auto index = fmindex_collection::ReverseFMIndex<OccTable, fmindex_collection::DenseCSA>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            INFO(i);
            INFO(std::get<1>(index.locate(i)));
            INFO(sa[i]);
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }

    }

    SECTION("compare to a directly created index") {
        auto bwt = std::vector<uint8_t>{'H', 'W', 'a', 'e', 'l', 'l', 'l', 't', 'o', ' ', '\0'};
        auto sa  = DenseVector{11, 6, 1, 7, 2, 8, 3, 9, 4, 5, 10};

        auto text  = std::vector<uint8_t>{'H', 'a', 'l', 'l', 'o', ' ', 'W', 'e', 'l', 't'};
        auto index = fmindex_collection::ReverseFMIndex<OccTable, fmindex_collection::DenseCSA>{std::vector<std::vector<uint8_t>>{text}, /*samplingRate*/1, /*threadNbr*/1};

        REQUIRE(bwt.size() == index.size());
        REQUIRE(sa.size() == index.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            INFO(i);
            INFO(sa[i]);
            INFO(std::get<1>(index.locate(i)));

            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }

    SECTION("compare to a directly created index but with sampling") {
        auto bwt = std::vector<uint8_t>{'H', '\0', 'W', 'a', 'e', 'l', 'l', 'l', 't', 'o', ' ', '\0'};
        auto sa  = DenseVector{0, 11, 6, 1, 7, 2, 8, 3, 9, 4, 5, 10};

        auto text  = std::vector<uint8_t>{'H', 'a', 'l', 'l', 'o', ' ', 'W', 'e', 'l', 't'};
        auto index = fmindex_collection::ReverseFMIndex<OccTable, fmindex_collection::DenseCSA>{std::vector<std::vector<uint8_t>>{text}, /*samplingRate*/2, /*threadNbr*/1};

        REQUIRE(bwt.size() == index.size());
        REQUIRE(sa.size() == index.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            INFO(i);
            INFO(sa[i]);
            if (index.occ.symbol(i) != 0) { // symbol 0 is not locatable in sampling above 1
                INFO(std::get<1>(index.locate(i)));
                CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
            }
        }
    }

}
