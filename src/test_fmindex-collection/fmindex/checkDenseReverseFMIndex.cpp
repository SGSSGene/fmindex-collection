// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include "../occtables/allTables.h"

#include <catch2/catch_all.hpp>
#include <fmindex-collection/fmindex/ReverseFMIndex.h>
#include <fmindex-collection/suffixarray/DenseCSA.h>

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
        auto csa = fmindex_collection::DenseCSA{sa, bitStack};
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

        auto csa = fmindex_collection::DenseCSA{sa2, bitStack};
        auto index = fmindex_collection::ReverseFMIndex<OccTable, fmindex_collection::DenseCSA>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            INFO(i);
            INFO(std::get<1>(index.locate(i)));
            INFO(sa[i]);
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
            auto res = index.single_locate_step(i);
            INFO(i);
            INFO(sa[i]);
            if (sa[i] == sa.size()-1 || i % 2 == 0) {
                REQUIRE(res);
                CHECK(*res == std::make_tuple(0, sa[i]));
            } else {
                CHECK(!res);
            }
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

        auto csa = fmindex_collection::DenseCSA{sa2, bitStack};
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

        auto csa = fmindex_collection::DenseCSA{sa2, bitStack};
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
        auto bwt = std::vector<uint8_t>{'H', 'W', 'a', 'e', 'l', 'l', 'l', 't', 'o', ' ', '\0'};
        auto sa  = DenseVector{0, 6, 1, 7, 2, 8, 3, 9, 4, 5, 10};

        auto text  = std::vector<uint8_t>{'H', 'a', 'l', 'l', 'o', ' ', 'W', 'e', 'l', 't'};
        auto index = fmindex_collection::ReverseFMIndex<OccTable, fmindex_collection::DenseCSA>{std::vector<std::vector<uint8_t>>{text}, /*samplingRate*/2, /*threadNbr*/1};

        REQUIRE(bwt.size() == index.size());
        REQUIRE(sa.size() == index.size());
        //!TODO why can't position 0 be located?
        for (size_t i{1}; i < sa.size(); ++i) {
            INFO(i);
            INFO(sa[i]);
            INFO(std::get<1>(index.locate(i)));
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }

}
