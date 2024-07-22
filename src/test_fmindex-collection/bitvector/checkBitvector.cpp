// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include <catch2/catch_all.hpp>
#include <fmindex-collection/utils.h>
#include <fstream>
#include <nanobench.h>
#include <reflect>
#include <thread>

#include "allBitVectors.h"

#if __has_include(<cxxabi.h>)
#include <cxxabi.h>
template <typename T>
static auto getName() {
    int     status;
    auto realname = abi::__cxa_demangle(typeid(T).name(), nullptr, nullptr, &status);
    auto str = std::string{realname};
    std::free(realname);
    return str;
}
#else
template <typename T>
static auto getName() {
    return std::string{reflect::type_name<T>()};
}
#endif

TEMPLATE_TEST_CASE("check bit vectors are working", "[BitVector]", ALLBITVECTORS) {
    using Vector = TestType;
    auto vector_name = getName<Vector>();
    INFO(vector_name);

    SECTION("short text") {
        auto text = std::vector<uint8_t>{0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1};

        auto vec = Vector{text};
        REQUIRE(vec.size() == text.size());

        SECTION("check that symbol() call works") {
            for (size_t i{0}; i < text.size(); ++i) {
                INFO(i);
                CHECK(vec.symbol(i) == bool(text.at(i)));
            }
        }

        SECTION("test complete vector on rank()") {
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
        }
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

        auto vec = Vector{text};
        REQUIRE(vec.size() == text.size());

        SECTION("check that symbol() call works") {
            for (size_t i{0}; i < text.size(); ++i) {
                INFO(i);
                CHECK(vec.symbol(i) == bool(text.at(i)));
            }
        }

        SECTION("test complete vector on rank()") {
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
        }
    }
    SECTION("serialization/deserialization") {
        auto input = std::vector<uint8_t>{0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1};
        SECTION("serialize") {
            auto ofs = std::ofstream{"temp_test_serialization"};
            auto vec = Vector{input};
            auto archive = cereal::BinaryOutputArchive{ofs};
            archive(vec);
        }
        SECTION("deserialize") {
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
    }

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        auto rng = ankerl::nanobench::Rng{};

        auto text = std::vector<uint8_t>{};
        #ifdef NDEBUG
        for (size_t i{0}; i<100'000'000; ++i) {
        #else
        for (size_t i{0}; i<100'000; ++i) {
        #endif
            text.push_back(rng.bounded(4) == 0);
        }

        bench.batch(text.size()).run(vector_name + "()", [&]() {
            auto vec = Vector{text};
            ankerl::nanobench::doNotOptimizeAway(vec);
        });
        auto vec = Vector{text};
        bench.batch(1).run(vector_name + "::symbol()", [&]() {
            auto v = vec.symbol(rng.bounded(text.size()));
            ankerl::nanobench::doNotOptimizeAway(v);
        });
        bench.run(vector_name + "::rank()", [&]() {
            auto v = vec.rank(rng.bounded(text.size()));
            ankerl::nanobench::doNotOptimizeAway(v);
        });
    }
}
