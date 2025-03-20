// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include <catch2/catch_all.hpp>
#include <cstdlib>
#include <fmindex-collection/utils.h>
#include <fstream>
#include <nanobench.h>

#include <pasta/bit_vector/bit_vector.hpp>
#include <pasta/bit_vector/support/rank.hpp>
#include <pasta/bit_vector/support/flat_rank.hpp>
#include <pasta/bit_vector/support/wide_rank.hpp>
#include <rank9.h>

#include "../BenchSize.h"
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
#include <reflect>
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

        SECTION("test complete vector on symbol()") {
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

        auto vec2 = Vector{};
        for (auto c : text) {
            vec2.push_back(c);
        }

        SECTION("check that symbol() call works") {
            for (size_t i{0}; i < text.size(); ++i) {
                INFO(i);
                CHECK(vec.symbol(i) == bool(text.at(i)));
                CHECK(vec2.symbol(i) == bool(text.at(i)));
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
    }

    SECTION("very longer text") {
        srand(0);
        auto text = std::vector<uint8_t>{};
        auto rank = std::vector<size_t>{0};
        for (size_t i{}; i < 65536ull*10; ++i) {
            text.push_back(rand()%2);
            rank.push_back(rank.back() + text.back());
        }

        auto vec = Vector{text};
        REQUIRE(vec.size() == text.size());

        SECTION("check that symbol() call works") {
            for (size_t i{0}; i < text.size(); ++i) {
                INFO(i);
                CHECK(vec.symbol(i) == bool(text.at(i)));
            }
        }

        SECTION("test complete vector on rank()") {
            for (size_t i{0}; i < text.size(); ++i) {
                INFO(i);
                CHECK(vec.rank(i) == rank.at(i));
            }
        }
    }

    SECTION("short text, push() back must have same result as c'tor") {
        srand(0);
        auto text = std::vector<uint8_t>{};
        auto rank = std::vector<size_t>{0};

        auto vec1 = Vector{};
        for (size_t i{}; i < 512; ++i) {
            text.push_back(rand()%2);
            rank.push_back(rank.back() + text.back());

            vec1.push_back(text.back());

            auto vec2 = Vector{text};
            SECTION("check vec2 and vec are the same") {
                for (size_t i{0}; i < text.size(); ++i) {
                    INFO(i);
                    CHECK(vec1.symbol(i) == bool(text.at(i)));
                    CHECK(vec2.symbol(i) == bool(text.at(i)));
                    CHECK(vec1.rank(i) == rank.at(i));
                    CHECK(vec2.rank(i) == rank.at(i));
                }
                CHECK(vec1.rank(text.size()) == rank.back());
                CHECK(vec2.rank(text.size()) == rank.back());


            }
        }
    }


    SECTION("short text, some inserted at creation, rest via push()") {
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
                auto v2 = vec.rank(vec.size());
                if (v1 != v2) {
                    auto vec2 = Vector{text};
                    CHECK(v1 == v2);

                }
                CHECK(v1 == v2);
            }
        }
        REQUIRE(vec.size() == text.size());

        SECTION("check that symbol() call works") {
            for (size_t i{0}; i < text.size(); ++i) {
                INFO(i);
                CHECK(vec.symbol(i) == bool(text.at(i)));
            }
        }

        SECTION("test complete vector on rank()") {
            for (size_t i{0}; i < text.size(); ++i) {
                INFO(i);
                CHECK(vec.rank(i) == rank.at(i));
            }
        }
    }


    SECTION("long text, some inserted at creation, rest via push()") {
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
                auto v2 = vec.rank(vec.size());
                CHECK(v1 == v2);
            }
        }
        REQUIRE(vec.size() == text.size());

        SECTION("check that symbol() call works") {
            for (size_t i{0}; i < text.size(); ++i) {
                INFO(i);
                CHECK(vec.symbol(i) == bool(text.at(i)));
            }
        }

        SECTION("test complete vector on rank()") {
            for (size_t i{0}; i < text.size(); ++i) {
                INFO(i);
                CHECK(vec.rank(i) == rank.at(i));
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
}

namespace {

// Wrapper for other bitvectors
struct FlatRank {
    pasta::BitVector bv;
    pasta::FlatRank<pasta::OptimizedFor::DONT_CARE, pasta::BitVector> rs{bv};
    FlatRank() = default;
    FlatRank(FlatRank&&) = default;
    FlatRank(std::vector<bool> const& input)
        : bv{[&]() {
            pasta::BitVector bv{input.size(), 0};
            for (size_t i = 0; i < bv.size(); ++i) {
                bv[i] = input[i];
            }
            return bv;
        }()}
    {
    }
    auto symbol(size_t i) const -> bool {
        return bv[i];
    }
    auto rank(size_t i) const -> size_t {
        return rs.rank1(i);
    }
    auto space_usage() const -> size_t {
        return bv.space_usage() + rs.space_usage();
    }

};

struct WideRank {
    pasta::BitVector bv;
    pasta::WideRank<pasta::OptimizedFor::DONT_CARE, pasta::BitVector> rs{bv};
    WideRank() = default;
    WideRank(WideRank&&) = default;
    WideRank(std::vector<bool> const& input)
        : bv{[&]() {
            pasta::BitVector bv{input.size(), 0};
            for (size_t i = 0; i < bv.size(); ++i) {
                bv[i] = input[i];
            }
            return bv;
        }()}
    {
    }
    auto symbol(size_t i) const -> bool {
        return bv[i];
    }
    auto rank(size_t i) const -> size_t {
        return rs.rank1(i);
    }
    auto space_usage() const -> size_t {
        return bv.space_usage() + rs.space_usage();
    }
};
struct Rank9 {
    std::vector<uint64_t> bitvector;
    rank9 bv;

    Rank9() = default;
    Rank9(Rank9&&) = default;
    Rank9(std::vector<bool> input)
        : bitvector{[&]() {
            auto bitvector = std::vector<uint64_t>{};
            bitvector.resize(input.size()/64 + 1);
            for (size_t i = 0; i < input.size(); ++i) {
                auto id = i / 64;
                auto offset = i % 64;
                bitvector[id] |= (input[i]<<offset);
            }
            return bitvector;
        }()}
        , bv{bitvector.data(), input.size()+1}
    {}

    auto symbol(size_t i) const -> bool {
        auto id = i / 64;
        auto offset = i % 64;
        return (bitvector[id] >> offset) & 1;
    }

    auto rank(size_t i) const -> size_t {
        //hack, since rank9 is not const correct
        return const_cast<rank9&>(bv).rank(i);
    }
    auto space_usage() const -> size_t {
        return 0;
//        bitvector.size() * 8 + rank9.
    }

};

template <typename ...T>
void call_with_templates(auto f) {
    (f.template operator()<T>(), ...);
}

auto generateText() -> std::vector<bool> const& {
    static auto text = []() -> std::vector<bool> {
        auto rng = ankerl::nanobench::Rng{};

        auto size = []() -> size_t {
            auto ptr = std::getenv("BITVECTORSIZE");
            if (ptr) {
                return std::stoull(ptr);
            }
            #ifdef NDEBUG
                return 10'000'000;
            #else
                return 100'000;
            #endif
        }();

        auto text = std::vector<bool>{};
        for (size_t i{0}; i<size; ++i) {
            text.push_back(rng.bounded(4) == 0);
        }
        return text;
    }();
    return text;
}
}

TEST_CASE("benchmark bit vectors ctor run times", "[BitVector][!benchmark][time][ctor][.]") {
    auto bench_ctor = ankerl::nanobench::Bench{};
    bench_ctor.title("c'tor()")
              .relative(true);

    auto& text = generateText();

    SECTION("benchmarking") {
        call_with_templates<
            ALLBITVECTORS,
            FlatRank,
            WideRank,
            Rank9>([&]<typename Vector>() {

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            bench_ctor.batch(text.size()).run(vector_name, [&]() {
                auto vec = Vector{text};
                ankerl::nanobench::doNotOptimizeAway(vec);
            });
        });
    }
}

TEST_CASE("benchmark bit vectors rank and symbol run times", "[BitVector][!benchmark][time][rank][symbol]") {

    auto& text = generateText();

    SECTION("benchmarking - symbol") {
        auto bench_symbol = ankerl::nanobench::Bench{};
        bench_symbol.title("symbol()")
                    .relative(true);

        bench_symbol.epochs(10);
        bench_symbol.minEpochTime(std::chrono::milliseconds{10});
        call_with_templates<
            ALLBITVECTORS,
            FlatRank,
            WideRank,
            Rank9>([&]<typename Vector>() {

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto rng = ankerl::nanobench::Rng{};

            auto vec = Vector{text};

            bench_symbol.run(vector_name, [&]() {
                auto v = vec.symbol(rng.bounded(text.size()));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        });
    }

    SECTION("benchmarking - rank") {
        auto bench_rank = ankerl::nanobench::Bench{};
        bench_rank.title("rank()")
                  .relative(true);

        bench_rank.epochs(20);
        bench_rank.minEpochTime(std::chrono::milliseconds{1});
        bench_rank.minEpochIterations(1'000'000);

        call_with_templates<
            ALLBITVECTORS,
            FlatRank,
            WideRank,
            Rank9>([&]<typename Vector>() {

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto rng = ankerl::nanobench::Rng{};

            auto vec = Vector{text};

            bench_rank.run(vector_name, [&]() {
                auto v = vec.rank(rng.bounded(text.size()));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        });
    }
}

TEST_CASE("benchmark bit vectors memory consumption", "[BitVector][!benchmark][size]") {
    BenchSize benchSize;

    SECTION("benchmarking") {
        call_with_templates<
            ALLBITVECTORS,
            FlatRank,
            WideRank,
            Rank9>([&]<typename Vector>() {

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto& text = generateText();

            auto vec = Vector{text};
            if constexpr (requires() { vec.space_usage(); }) {
                auto s = vec.space_usage();
                benchSize.addEntry({
                    .name = vector_name,
                    .size = s,
                    .text_size = text.size(),
                    .bits_per_char = (s*8)/double(text.size())
                });
            } else {
                auto ofs     = std::stringstream{};
                auto archive = cereal::BinaryOutputArchive{ofs};
                archive(vec);
                auto s = ofs.str().size();
                benchSize.addEntry({
                    .name = vector_name,
                    .size = s,
                    .text_size = text.size(),
                    .bits_per_char = (s*8)/double(text.size())
                });
            }
        });
    }
}
