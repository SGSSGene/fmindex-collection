// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include <catch2/catch_all.hpp>
#include <cstdlib>
#include <fmindex-collection/utils.h>
#include <fstream>

#include "../string/utils.h"
#include "allBitVectors.h"

TEST_CASE("check bit vectors are working", "[bitvector]") {
    SECTION("short text") {
        auto text = std::vector<uint8_t>{0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1};

        call_with_templates<
            ALLBITVECTORS>([&]<typename Vector>() {
            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto vec = Vector{text};
            REQUIRE(vec.size() == text.size());

            // check that symbol() call works
            for (size_t i{0}; i < text.size(); ++i) {
                INFO(i);
                CHECK(vec.symbol(i) == bool(text.at(i)));
            }

            // test complete vector on symbol()
            CHECK(vec.symbol( 0) == 0);
            CHECK(vec.symbol( 1) == 1);
            CHECK(vec.symbol( 2) == 1);
            CHECK(vec.symbol( 3) == 0);
            CHECK(vec.symbol( 4) == 0);
            CHECK(vec.symbol( 5) == 1);
            CHECK(vec.symbol( 6) == 0);
            CHECK(vec.symbol( 7) == 1);
            CHECK(vec.symbol( 8) == 1);
            CHECK(vec.symbol( 9) == 1);
            CHECK(vec.symbol(10) == 0);
            CHECK(vec.symbol(11) == 0);
            CHECK(vec.symbol(12) == 0);
            CHECK(vec.symbol(13) == 1);

            // test complete vector on rank()
            CHECK(vec.rank( 0) == 0);
            CHECK(vec.rank( 1) == 0);
            CHECK(vec.rank( 2) == 1);
            CHECK(vec.rank( 3) == 2);
            CHECK(vec.rank( 4) == 2);
            CHECK(vec.rank( 5) == 2);
            CHECK(vec.rank( 6) == 3);
            CHECK(vec.rank( 7) == 3);
            CHECK(vec.rank( 8) == 4);
            CHECK(vec.rank( 9) == 5);
            CHECK(vec.rank(10) == 6);
            CHECK(vec.rank(11) == 6);
            CHECK(vec.rank(12) == 6);
            CHECK(vec.rank(13) == 6);
            CHECK(vec.rank(14) == 7);
        });
    }
    SECTION("longer text") {
        auto text = std::vector<uint8_t>{0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
                                         0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,

                                        };

        call_with_templates<
            ALLBITVECTORS>([&]<typename Vector>() {
            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto vec = Vector{text};
            REQUIRE(vec.size() == text.size());

            // check that symbol() call works
            for (size_t i{0}; i < text.size(); ++i) {
                INFO(i);
                CHECK(vec.symbol(i) == bool(text.at(i)));
            }

            // test complete vector on rank()
            for (size_t i{0}; i < 512; i += 16) {
                CHECK(vec.rank(i +  0) == 0 + i/16*8);
                CHECK(vec.rank(i +  1) == 0 + i/16*8);
                CHECK(vec.rank(i +  2) == 1 + i/16*8);
                CHECK(vec.rank(i +  3) == 2 + i/16*8);
                CHECK(vec.rank(i +  4) == 2 + i/16*8);
                CHECK(vec.rank(i +  5) == 2 + i/16*8);
                CHECK(vec.rank(i +  6) == 3 + i/16*8);
                CHECK(vec.rank(i +  7) == 3 + i/16*8);
                CHECK(vec.rank(i +  8) == 4 + i/16*8);
                CHECK(vec.rank(i +  9) == 5 + i/16*8);
                CHECK(vec.rank(i + 10) == 6 + i/16*8);
                CHECK(vec.rank(i + 11) == 6 + i/16*8);
                CHECK(vec.rank(i + 12) == 6 + i/16*8);
                CHECK(vec.rank(i + 13) == 6 + i/16*8);
                CHECK(vec.rank(i + 14) == 7 + i/16*8);
                CHECK(vec.rank(i + 15) == 8 + i/16*8);
                CHECK(vec.rank(i + 16) == 8 + i/16*8);
            }

            // check push back is working (if push_back available)
            if constexpr (requires { Vector{}.push_back(uint8_t{0}); }) {
                auto vec2 = Vector{};
                for (auto c : text) {
                    vec2.push_back(c);
                }
                REQUIRE(vec2.size() == text.size());

                // check that symbol() call works
                for (size_t i{0}; i < text.size(); ++i) {
                    INFO(i);
                    CHECK(vec2.symbol(i) == bool(text.at(i)));
                }
                // test complete vector on rank()
                for (size_t i{0}; i < 512; i += 16) {
                    CHECK(vec2.rank(i +  0) == 0 + i/16*8);
                    CHECK(vec2.rank(i +  1) == 0 + i/16*8);
                    CHECK(vec2.rank(i +  2) == 1 + i/16*8);
                    CHECK(vec2.rank(i +  3) == 2 + i/16*8);
                    CHECK(vec2.rank(i +  4) == 2 + i/16*8);
                    CHECK(vec2.rank(i +  5) == 2 + i/16*8);
                    CHECK(vec2.rank(i +  6) == 3 + i/16*8);
                    CHECK(vec2.rank(i +  7) == 3 + i/16*8);
                    CHECK(vec2.rank(i +  8) == 4 + i/16*8);
                    CHECK(vec2.rank(i +  9) == 5 + i/16*8);
                    CHECK(vec2.rank(i + 10) == 6 + i/16*8);
                    CHECK(vec2.rank(i + 11) == 6 + i/16*8);
                    CHECK(vec2.rank(i + 12) == 6 + i/16*8);
                    CHECK(vec2.rank(i + 13) == 6 + i/16*8);
                    CHECK(vec2.rank(i + 14) == 7 + i/16*8);
                    CHECK(vec2.rank(i + 15) == 8 + i/16*8);
                    CHECK(vec2.rank(i + 16) == 8 + i/16*8);
                }
            }
        });
    }

    SECTION("very longer text") {
        call_with_templates<
            ALLBITVECTORS>([&]<typename Vector>() {
            auto vector_name = getName<Vector>();
            INFO(vector_name);

            srand(0);
            auto text = std::vector<uint8_t>{};
            auto rank = std::vector<size_t>{0};
            for (size_t i{}; i < 65536ull*10; ++i) {
                text.push_back(rand()%2);
                rank.push_back(rank.back() + text.back());
            }

            auto vec = Vector{text};
            REQUIRE(vec.size() == text.size());

            // check that symbol() call works
            for (size_t i{0}; i < text.size(); ++i) {
                INFO(i);
                CHECK(vec.symbol(i) == bool(text.at(i)));
            }

            // test complete vector on rank()
            for (size_t i{0}; i < text.size(); ++i) {
                INFO(i);
                CHECK(vec.rank(i) == rank.at(i));
            }
        });
    }

    SECTION("short text, push() back must have same result as c'tor") {
        call_with_templates<
            ALLBITVECTORS>([&]<typename Vector>() {
            auto vector_name = getName<Vector>();
            INFO(vector_name);

            // check only if push_back is available
            if constexpr (requires { Vector{}.push_back(uint8_t{0}); }) {
                srand(0);
                auto text = std::vector<uint8_t>{};
                auto rank = std::vector<size_t>{0};

                auto vec1 = Vector{};
                for (size_t i{}; i < 512; ++i) {
                    text.push_back(rand()%2);
                    rank.push_back(rank.back() + text.back());

                    vec1.push_back(text.back());

                    auto vec2 = Vector{text};
                    auto vec3 = vec1;
                    // check vec2 and vec are the same
                    for (size_t i{0}; i < text.size(); ++i) {
                        INFO(i);
                        CHECK(vec3.symbol(i) == bool(text.at(i)));
                        CHECK(vec2.symbol(i) == bool(text.at(i)));
                        CHECK(vec3.rank(i) == rank.at(i));
                        CHECK(vec2.rank(i) == rank.at(i));
                    }
                    CHECK(vec3.rank(text.size()) == rank.back());
                    CHECK(vec2.rank(text.size()) == rank.back());
                }
            }
        });
    }

    SECTION("short text, some inserted at creation, rest via push()") {
        call_with_templates<
            ALLBITVECTORS>([&]<typename Vector>() {
            auto vector_name = getName<Vector>();
            INFO(vector_name);
            // check only if push_back is available
            if constexpr (requires { Vector{}.push_back(uint8_t{0}); }) {
                srand(0);
                auto text = std::vector<uint8_t>{};
                auto rank = std::vector<size_t>{0};
                for (size_t i{}; i < 64; ++i) {
                    text.push_back(rand()%2);
                    rank.push_back(rank.back() + text.back());
                }

                auto vec = Vector{text};
                REQUIRE(vec.size() == text.size());

                for (size_t i{}; i < 64; ++i) {
                    text.push_back(rand()%2);
                    rank.push_back(rank.back() + text.back());
                    vec.push_back(text.back());
                    {
                        auto v1 = rank.back();
                        auto vec2 = vec;
                        auto v2 = vec2.rank(vec.size());
                        CHECK(v1 == v2);
                    }
                }
                REQUIRE(vec.size() == text.size());

                // check that symbol() call works
                for (size_t i{0}; i < text.size(); ++i) {
                    INFO(i);
                    CHECK(vec.symbol(i) == bool(text.at(i)));
                }

                // test complete vector on rank()
                for (size_t i{0}; i < text.size(); ++i) {
                    INFO(i);
                    CHECK(vec.rank(i) == rank.at(i));
                }
            }
        });
    }

    SECTION("long text, some inserted at creation, rest via push()") {
        call_with_templates<
            ALLBITVECTORS>([&]<typename Vector>() {
            auto vector_name = getName<Vector>();
            INFO(vector_name);

            // check only if push_back is available
            if constexpr (requires { Vector{}.push_back(uint8_t{0}); }) {
                srand(0);
                auto text = std::vector<uint8_t>{};
                auto rank = std::vector<size_t>{0};
                for (size_t i{}; i < (size_t{1ul}<<15); ++i) {
                    text.push_back(rand()%2);
                    rank.push_back(rank.back() + text.back());
                }

                auto vec = Vector{text};
                REQUIRE(vec.size() == text.size());

                for (size_t i{}; i < (size_t{1ul}<<15); ++i) {
                    text.push_back(rand()%2);
                    rank.push_back(rank.back() + text.back());
                    vec.push_back(text.back());
                    {
                        auto v1 = rank.back();
                        auto vec3 = vec;
                        auto v2 = vec3.rank(vec.size());
                        CHECK(v1 == v2);
                    }
                }
                REQUIRE(vec.size() == text.size());

                // check that symbol() call works
                for (size_t i{0}; i < text.size(); ++i) {
                    INFO(i);
                    CHECK(vec.symbol(i) == bool(text.at(i)));
                }

                // test complete vector on rank()
                for (size_t i{0}; i < text.size(); ++i) {
                    INFO(i);
                    CHECK(vec.rank(i) == rank.at(i));
                }
            }
        });
    }


    SECTION("serialization/deserialization") {
        call_with_templates<
            ALLBITVECTORS>([&]<typename Vector>() {
            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto input = std::vector<uint8_t>{0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1};
            // serialize
            {
                auto ofs = std::ofstream{"temp_test_serialization"};
                auto vec = Vector{input};
                auto archive = cereal::BinaryOutputArchive{ofs};
                archive(vec);
            }
            // deserialize
            {
                auto ifs = std::ifstream{"temp_test_serialization"};
                auto vec = Vector{};
                auto archive = cereal::BinaryInputArchive{ifs};
                archive(vec);
                REQUIRE(input.size() == vec.size());
                size_t count{};
                for (size_t i{0}; i != input.size(); ++i) {
                    CHECK((bool)input[i] == vec.symbol(i));
                    CHECK(count == vec.rank(i));
                    count += input[i];
                }
                CHECK(count == vec.rank(input.size()));
            }
        });
    }
}
