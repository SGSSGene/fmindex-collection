// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include "../occtables/allTables.h"

#include <catch2/catch_all.hpp>
#include <fmindex-collection/fmindex/BiFMIndex.h>
#include <fstream>

TEMPLATE_TEST_CASE("checking bidirectional fm index", "[BiFMIndex]", ALLTABLES) {
    using OccTable = TestType;

    auto bwt    = std::vector<uint8_t>{'t', '\0', 'o', '\0', ' ', 'H', 'W', 'a', 'l', 'e', 'l', 'l'};
    auto bwtRev = std::vector<uint8_t>{'H', '\0', 'W', 'a', 'e', 'l', 'l', 'l', 't', 'o', ' ', '\0'};
    auto sa     = std::vector<uint64_t>{ 10, 11, 5, 0,  6,  1,  7,  2,  3,  8,  4,  9 };

    SECTION("full sa") {
        auto bitStack = fmindex_collection::BitStack{};
        for (size_t i{0}; i < sa.size(); ++i) {
            bitStack.push(true);
        }
        auto csa = fmindex_collection::CSA{sa, bitStack, /*.threadNbr=*/63, /*.seqCount=*/1};
        auto index = fmindex_collection::BiFMIndex<OccTable>{bwt, bwtRev, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
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

        auto csa = fmindex_collection::CSA{sa2, bitStack, /*.threadNbr=*/63, /*.seqCount=*/1};
        auto index = fmindex_collection::BiFMIndex<OccTable>{bwt, bwtRev, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
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
        auto bitStack = fmindex_collection::BitStack{};
        auto sa2 = std::vector<uint64_t>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            auto add = bool{i % 2 == 1};
            bitStack.push(add);
            if (add) {
                sa2.push_back(sa[i]);
            }
        }

        auto csa = fmindex_collection::CSA{sa2, bitStack, /*.threadNbr=*/63, /*.seqCount=*/1};
        auto index = fmindex_collection::BiFMIndex<OccTable>{bwt, bwtRev, std::move(csa)};

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

        auto csa = fmindex_collection::CSA{sa2, bitStack, /*.threadNbr=*/63, /*.seqCount=*/1};
        auto index = fmindex_collection::BiFMIndex<OccTable>{bwt, bwtRev, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }


    SECTION("serialization/deserialization") {
        SECTION("serialize") {
            auto ofs = std::ofstream{"temp_test_serialization", std::ios::binary};

            auto bitStack = fmindex_collection::BitStack{};
            for (size_t i{0}; i < sa.size(); ++i) {
                bitStack.push(true);
            }
            auto csa = fmindex_collection::CSA{sa, bitStack, /*.threadNbr=*/63, /*.seqCount=*/1};
            auto index = fmindex_collection::BiFMIndex<OccTable>{bwt, bwtRev, std::move(csa)};
            auto archive = cereal::BinaryOutputArchive{ofs};
            archive(index);
        }
        SECTION("deserialize") {
            auto ifs = std::ifstream{"temp_test_serialization", std::ios::binary};

            auto index = fmindex_collection::BiFMIndex<OccTable>{};
            auto archive = cereal::BinaryInputArchive{ifs};
            archive(index);

            REQUIRE(index.size() == bwt.size());
            for (size_t i{0}; i < sa.size(); ++i) {
                CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
            }
        }
    }

}

TEMPLATE_TEST_CASE("checking bidirectional fm index on longer text (more than 256 chars)", "[BiFMIndex]", ALLTABLES) {
    using OccTable = TestType;

    auto bwt    = std::vector<uint8_t>{'0', '4', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '4', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '2', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '4', '4', '4', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '0', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '3', '4', '3', '3', '3', '4', '4', '4', '\0', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8' };
    auto bwtRev = std::vector<uint8_t>{'4', '\0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '2', '1', '1', '1', '1', '1', '1', '1', '1', '3', '1', '1', '1', '1', '1', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '1', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '1', '3', '3', '3', '3', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '4', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '5', '2', '5', '5', '5', '5', '5', '5', '2', '2', '0', '4', '4', '4', '1', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '8', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' };
    auto sa     = std::vector<uint64_t>{329, 328, 314, 304, 294, 284, 274, 264, 254, 244, 234, 224, 214, 204, 194, 184, 174, 164, 75, 85, 95, 105, 17, 115, 27, 125, 37, 135, 47, 145, 57, 155, 67, 5, 6, 315, 305, 295, 285, 275, 265, 255, 245, 235, 225, 215, 205, 195, 185, 175, 165, 76, 86, 96, 8, 106, 18, 116, 28, 126, 38, 136, 48, 146, 58, 7, 325, 322, 319, 316, 306, 296, 286, 276, 266, 256, 246, 236, 226, 216, 206, 196, 186, 176, 166, 156, 77, 87, 97, 9, 107, 19, 117, 29, 127, 39, 137, 49, 147, 59, 326, 323, 320, 317, 307, 297, 287, 277, 267, 257, 247, 237, 227, 217, 207, 197, 187, 177, 167, 157, 68, 78, 88, 98, 10, 108, 20, 118, 30, 128, 40, 138, 50, 148, 60, 327, 4, 324, 321, 318, 3, 2, 1, 0, 308, 298, 288, 278, 268, 258, 248, 238, 228, 218, 208, 198, 188, 178, 168, 158, 69, 79, 89, 99, 11, 109, 21, 119, 31, 129, 41, 139, 51, 149, 61, 309, 299, 289, 279, 269, 259, 249, 239, 229, 219, 209, 199, 189, 179, 169, 159, 70, 80, 90, 100, 12, 110, 22, 120, 32, 130, 42, 140, 52, 150, 62, 310, 300, 290, 280, 270, 260, 250, 240, 230, 220, 210, 200, 190, 180, 170, 160, 71, 81, 91, 101, 13, 111, 23, 121, 33, 131, 43, 141, 53, 151, 63, 311, 301, 291, 281, 271, 261, 251, 241, 231, 221, 211, 201, 191, 181, 171, 161, 72, 82, 92, 102, 14, 112, 24, 122, 34, 132, 44, 142, 54, 152, 64, 312, 302, 292, 282, 272, 262, 252, 242, 232, 222, 212, 202, 192, 182, 172, 162, 73, 83, 93, 103, 15, 113, 25, 123, 35, 133, 45, 143, 55, 153, 65, 313, 303, 293, 283, 273, 263, 253, 243, 233, 223, 213, 203, 193, 183, 173, 163, 74, 84, 94, 104, 16, 114, 26, 124, 36, 134, 46, 144, 56, 154, 66, };

    SECTION("full sa") {
        auto bitStack = fmindex_collection::BitStack{};
        for (size_t i{0}; i < sa.size(); ++i) {
            bitStack.push(true);
        }
        auto csa = fmindex_collection::CSA{sa, bitStack, /*.threadNbr=*/63, /*.seqCount=*/1};
        auto index = fmindex_collection::BiFMIndex<OccTable>{bwt, bwtRev, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
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

        auto csa = fmindex_collection::CSA{sa2, bitStack, /*.threadNbr=*/63, /*.seqCount=*/1};
        auto index = fmindex_collection::BiFMIndex<OccTable>{bwt, bwtRev, std::move(csa)};

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

        auto csa = fmindex_collection::CSA{sa2, bitStack, /*.threadNbr=*/63, /*.seqCount=*/1};
        auto index = fmindex_collection::BiFMIndex<OccTable>{bwt, bwtRev, std::move(csa)};

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

        auto csa = fmindex_collection::CSA{sa2, bitStack, /*.threadNbr=*/63, /*.seqCount=*/1};
        auto index = fmindex_collection::BiFMIndex<OccTable>{bwt, bwtRev, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }
}
