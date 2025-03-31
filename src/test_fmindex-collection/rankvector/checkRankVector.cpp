// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include <catch2/catch_all.hpp>
#include <fmindex-collection/utils.h>
#include <fstream>
#include <nanobench.h>

#include "../BenchSize.h"
#include "allRankVectors.h"

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

namespace {

template <typename ...T>
void call_with_templates(auto f) {
    (f.template operator()<T>(), ...);
}

template <size_t min=1, size_t range=4>
auto generateText() -> std::vector<uint8_t> const& {
    static auto text = []() -> std::vector<uint8_t> {
        auto rng = ankerl::nanobench::Rng{};

        // generates string with values between 1-4
        auto size = []() -> size_t {
            auto ptr = std::getenv("VECTORSIZE");
            if (ptr) {
                return std::stoull(ptr);
            }
            #ifdef NDEBUG
                return 1'000'000;
            #else
                return 1'000;
            #endif
        }();

        auto text = std::vector<uint8_t>{};
        for (size_t i{0}; i<size; ++i) {
            text.push_back(rng.bounded(range) + min);
        }
        return text;
    }();
    return text;
}


struct Bench : ankerl::nanobench::Bench {
    std::stringstream output{};

    Bench(std::string title) {
        this->title(title)
            .relative(true)
            .output(&output)
        ;
    }
    ~Bench() {
        if (output.str().size() > 0) {
            std::cout << output.str() << '\n';
        }
    }
};

struct Benchs {
    Bench bench_rank{"rank()"};
    Bench bench_prefix_rank{"prefix_rank()"};
    Bench bench_all_ranks{"all_ranks()"};
    Bench bench_all_prefix_ranks{"all_prefix_ranks()"};
    Bench bench_symbol{"symbol()"};
    Bench bench_ctor{"c'tor"};
};
}


#if 0
TEMPLATE_TEST_CASE("check if rank on the symbol vectors is working", "[RankVector][256]", ALLRANKVECTORS(256)) {
    using Vector = TestType;
    auto vector_name = getName<Vector>();
    INFO(vector_name);

    auto text = std::vector<uint8_t>{'H', 'a', 'l', 'l', 'o', ' ', 'W', 'e', 'l', 't'};

    auto vec = Vector{std::span{text}};
    REQUIRE(vec.size() == text.size());

    SECTION("check that symbol() call works") {
        for (size_t i{0}; i < text.size(); ++i) {
            INFO(i);
            CHECK(vec.symbol(i) == text.at(i));
        }
    }
    CHECK(Vector::Sigma >= 128);

    SECTION("test complete vector on symbol()") {
        CHECK(vec.symbol( 0) == 'H');
        CHECK(vec.symbol( 1) == 'a');
        CHECK(vec.symbol( 2) == 'l');
        CHECK(vec.symbol( 3) == 'l');
        CHECK(vec.symbol( 4) == 'o');
        CHECK(vec.symbol( 5) == ' ');
        CHECK(vec.symbol( 6) == 'W');
        CHECK(vec.symbol( 7) == 'e');
        CHECK(vec.symbol( 8) == 'l');
        CHECK(vec.symbol( 9) == 't');
    }

    SECTION("test complete vector on rank()") {
        CHECK(vec.rank( 0, ' ') == 0);
        CHECK(vec.rank( 1, ' ') == 0);
        CHECK(vec.rank( 2, ' ') == 0);
        CHECK(vec.rank( 3, ' ') == 0);
        CHECK(vec.rank( 4, ' ') == 0);
        CHECK(vec.rank( 5, ' ') == 0);
        CHECK(vec.rank( 6, ' ') == 1);
        CHECK(vec.rank( 7, ' ') == 1);
        CHECK(vec.rank( 8, ' ') == 1);
        CHECK(vec.rank( 9, ' ') == 1);
        CHECK(vec.rank(10, ' ') == 1);

        CHECK(vec.rank( 0, 'H') == 0);
        CHECK(vec.rank( 1, 'H') == 1);
        CHECK(vec.rank( 2, 'H') == 1);
        CHECK(vec.rank( 3, 'H') == 1);
        CHECK(vec.rank( 4, 'H') == 1);
        CHECK(vec.rank( 5, 'H') == 1);
        CHECK(vec.rank( 6, 'H') == 1);
        CHECK(vec.rank( 7, 'H') == 1);
        CHECK(vec.rank( 8, 'H') == 1);
        CHECK(vec.rank( 9, 'H') == 1);
        CHECK(vec.rank(10, 'H') == 1);

        CHECK(vec.rank( 0, 'W') == 0);
        CHECK(vec.rank( 1, 'W') == 0);
        CHECK(vec.rank( 2, 'W') == 0);
        CHECK(vec.rank( 3, 'W') == 0);
        CHECK(vec.rank( 4, 'W') == 0);
        CHECK(vec.rank( 5, 'W') == 0);
        CHECK(vec.rank( 6, 'W') == 0);
        CHECK(vec.rank( 7, 'W') == 1);
        CHECK(vec.rank( 8, 'W') == 1);
        CHECK(vec.rank( 9, 'W') == 1);
        CHECK(vec.rank(10, 'W') == 1);

        CHECK(vec.rank( 0, 'a') == 0);
        CHECK(vec.rank( 1, 'a') == 0);
        CHECK(vec.rank( 2, 'a') == 1);
        CHECK(vec.rank( 3, 'a') == 1);
        CHECK(vec.rank( 4, 'a') == 1);
        CHECK(vec.rank( 5, 'a') == 1);
        CHECK(vec.rank( 6, 'a') == 1);
        CHECK(vec.rank( 7, 'a') == 1);
        CHECK(vec.rank( 8, 'a') == 1);
        CHECK(vec.rank( 9, 'a') == 1);
        CHECK(vec.rank(10, 'a') == 1);

        CHECK(vec.rank( 0, 'e') == 0);
        CHECK(vec.rank( 1, 'e') == 0);
        CHECK(vec.rank( 2, 'e') == 0);
        CHECK(vec.rank( 3, 'e') == 0);
        CHECK(vec.rank( 4, 'e') == 0);
        CHECK(vec.rank( 5, 'e') == 0);
        CHECK(vec.rank( 6, 'e') == 0);
        CHECK(vec.rank( 7, 'e') == 0);
        CHECK(vec.rank( 8, 'e') == 1);
        CHECK(vec.rank( 9, 'e') == 1);
        CHECK(vec.rank(10, 'e') == 1);

        CHECK(vec.rank( 0, 'l') == 0);
        CHECK(vec.rank( 1, 'l') == 0);
        CHECK(vec.rank( 2, 'l') == 0);
        CHECK(vec.rank( 3, 'l') == 1);
        CHECK(vec.rank( 4, 'l') == 2);
        CHECK(vec.rank( 5, 'l') == 2);
        CHECK(vec.rank( 6, 'l') == 2);
        CHECK(vec.rank( 7, 'l') == 2);
        CHECK(vec.rank( 8, 'l') == 2);
        CHECK(vec.rank( 9, 'l') == 3);
        CHECK(vec.rank(10, 'l') == 3);

        CHECK(vec.rank( 0, 'o') == 0);
        CHECK(vec.rank( 1, 'o') == 0);
        CHECK(vec.rank( 2, 'o') == 0);
        CHECK(vec.rank( 3, 'o') == 0);
        CHECK(vec.rank( 4, 'o') == 0);
        CHECK(vec.rank( 5, 'o') == 1);
        CHECK(vec.rank( 6, 'o') == 1);
        CHECK(vec.rank( 7, 'o') == 1);
        CHECK(vec.rank( 8, 'o') == 1);
        CHECK(vec.rank( 9, 'o') == 1);
        CHECK(vec.rank(10, 'o') == 1);

        CHECK(vec.rank( 0, 't') == 0);
        CHECK(vec.rank( 1, 't') == 0);
        CHECK(vec.rank( 2, 't') == 0);
        CHECK(vec.rank( 3, 't') == 0);
        CHECK(vec.rank( 4, 't') == 0);
        CHECK(vec.rank( 5, 't') == 0);
        CHECK(vec.rank( 6, 't') == 0);
        CHECK(vec.rank( 7, 't') == 0);
        CHECK(vec.rank( 8, 't') == 0);
        CHECK(vec.rank( 9, 't') == 0);
        CHECK(vec.rank(10, 't') == 1);

        // check all other characters have rank 0
        auto ignore = std::unordered_set<size_t>{' ', 'H', 'W', 'a', 'e', 'l', 'o', 't'};
        for (size_t s{0}; s < 256; ++s) {
            if (ignore.contains(s)) continue;
            for (size_t i{0}; i < 11; ++i) {
                CHECK(vec.rank(i, s) == 0);
            }
        }
    }

    SECTION("test complete vec 'H' for prefix_rank()") {
        CHECK(vec.prefix_rank( 0, ' ') == 0);
        CHECK(vec.prefix_rank( 1, ' ') == 0);
        CHECK(vec.prefix_rank( 2, ' ') == 0);
        CHECK(vec.prefix_rank( 3, ' ') == 0);
        CHECK(vec.prefix_rank( 4, ' ') == 0);
        CHECK(vec.prefix_rank( 5, ' ') == 0);
        size_t a{};
        for (size_t s{0}; s < ' '; ++s) {
            a += vec.rank(6, s);
            INFO(a);
            CHECK(a == 0);
        }
        CHECK(vec.rank(6, ' ' ) == 1);
        CHECK(vec.prefix_rank( 6, ' ') == 1);
        CHECK(vec.prefix_rank( 7, ' ') == 1);
        CHECK(vec.prefix_rank( 8, ' ') == 1);
        CHECK(vec.prefix_rank( 9, ' ') == 1);
        CHECK(vec.prefix_rank(10, ' ') == 1);

        CHECK(vec.prefix_rank( 0, 'H') == 0);
        CHECK(vec.prefix_rank( 1, 'H') == 1);
        CHECK(vec.prefix_rank( 2, 'H') == 1);
        CHECK(vec.prefix_rank( 3, 'H') == 1);
        CHECK(vec.prefix_rank( 4, 'H') == 1);
        CHECK(vec.prefix_rank( 5, 'H') == 1);
        CHECK(vec.prefix_rank( 6, 'H') == 2);
        CHECK(vec.prefix_rank( 7, 'H') == 2);
        CHECK(vec.prefix_rank( 8, 'H') == 2);
        CHECK(vec.prefix_rank( 9, 'H') == 2);
        CHECK(vec.prefix_rank(10, 'H') == 2);

        CHECK(vec.prefix_rank( 0, 'W') == 0);
        CHECK(vec.prefix_rank( 1, 'W') == 1);
        CHECK(vec.prefix_rank( 2, 'W') == 1);
        CHECK(vec.prefix_rank( 3, 'W') == 1);
        CHECK(vec.prefix_rank( 4, 'W') == 1);
        CHECK(vec.prefix_rank( 5, 'W') == 1);
        CHECK(vec.prefix_rank( 6, 'W') == 2);
        CHECK(vec.prefix_rank( 7, 'W') == 3);
        CHECK(vec.prefix_rank( 8, 'W') == 3);
        CHECK(vec.prefix_rank( 9, 'W') == 3);
        CHECK(vec.prefix_rank(10, 'W') == 3);

        CHECK(vec.prefix_rank( 0, 'a') == 0);
        CHECK(vec.prefix_rank( 1, 'a') == 1);
        CHECK(vec.prefix_rank( 2, 'a') == 2);
        CHECK(vec.prefix_rank( 3, 'a') == 2);
        CHECK(vec.prefix_rank( 4, 'a') == 2);
        CHECK(vec.prefix_rank( 5, 'a') == 2);
        CHECK(vec.prefix_rank( 6, 'a') == 3);
        CHECK(vec.prefix_rank( 7, 'a') == 4);
        CHECK(vec.prefix_rank( 8, 'a') == 4);
        CHECK(vec.prefix_rank( 9, 'a') == 4);
        CHECK(vec.prefix_rank(10, 'a') == 4);

        CHECK(vec.prefix_rank( 0, 'e') == 0);
        CHECK(vec.prefix_rank( 1, 'e') == 1);
        CHECK(vec.prefix_rank( 2, 'e') == 2);
        CHECK(vec.prefix_rank( 3, 'e') == 2);
        CHECK(vec.prefix_rank( 4, 'e') == 2);
        CHECK(vec.prefix_rank( 5, 'e') == 2);
        CHECK(vec.prefix_rank( 6, 'e') == 3);
        CHECK(vec.prefix_rank( 7, 'e') == 4);
        CHECK(vec.prefix_rank( 8, 'e') == 5);
        CHECK(vec.prefix_rank( 9, 'e') == 5);
        CHECK(vec.prefix_rank(10, 'e') == 5);

        CHECK(vec.prefix_rank( 0, 'l') == 0);
        CHECK(vec.prefix_rank( 1, 'l') == 1);
        CHECK(vec.prefix_rank( 2, 'l') == 2);
        CHECK(vec.prefix_rank( 3, 'l') == 3);
        CHECK(vec.prefix_rank( 4, 'l') == 4);
        CHECK(vec.prefix_rank( 5, 'l') == 4);
        CHECK(vec.prefix_rank( 6, 'l') == 5);
        CHECK(vec.prefix_rank( 7, 'l') == 6);
        CHECK(vec.prefix_rank( 8, 'l') == 7);
        CHECK(vec.prefix_rank( 9, 'l') == 8);
        CHECK(vec.prefix_rank(10, 'l') == 8);

        CHECK(vec.prefix_rank( 0, 'o') == 0);
        CHECK(vec.prefix_rank( 1, 'o') == 1);
        CHECK(vec.prefix_rank( 2, 'o') == 2);
        CHECK(vec.prefix_rank( 3, 'o') == 3);
        CHECK(vec.prefix_rank( 4, 'o') == 4);
        CHECK(vec.prefix_rank( 5, 'o') == 5);
        CHECK(vec.prefix_rank( 6, 'o') == 6);
        CHECK(vec.prefix_rank( 7, 'o') == 7);
        CHECK(vec.prefix_rank( 8, 'o') == 8);
        CHECK(vec.prefix_rank( 9, 'o') == 9);
        CHECK(vec.prefix_rank(10, 'o') == 9);

        CHECK(vec.prefix_rank( 0, 't') ==  0);
        CHECK(vec.prefix_rank( 1, 't') ==  1);
        CHECK(vec.prefix_rank( 2, 't') ==  2);
        CHECK(vec.prefix_rank( 3, 't') ==  3);
        CHECK(vec.prefix_rank( 4, 't') ==  4);
        CHECK(vec.prefix_rank( 5, 't') ==  5);
        CHECK(vec.prefix_rank( 6, 't') ==  6);
        CHECK(vec.prefix_rank( 7, 't') ==  7);
        CHECK(vec.prefix_rank( 8, 't') ==  8);
        CHECK(vec.prefix_rank( 9, 't') ==  9);
        CHECK(vec.prefix_rank(10, 't') == 10);
    }

    SECTION("check all_ranks() is equal to prefix_rank() and rank()") {
        for (size_t idx{0}; idx < vec.size(); ++idx) {
            auto [rank, prefix] = vec.all_ranks_and_prefix_ranks(idx);
            auto rank2 = vec.all_ranks(idx);
            for (size_t symb{1}; symb < Vector::Sigma; ++symb) {
                INFO(idx);
                INFO(symb);
                CHECK(rank[symb] == vec.rank(idx, symb));
                CHECK(rank2[symb] == vec.rank(idx, symb));
                CHECK(prefix[symb] == vec.prefix_rank(idx, symb));
            }
        }
    }
}

TEMPLATE_TEST_CASE("check symbol vectors construction on text longer than 256 characters", "[RankVector]", ALLRANKVECTORS(256)) {
    using Vector = TestType;
    auto vector_name = getName<Vector>();
    INFO(vector_name);

    auto text = std::vector<uint8_t>{'H', 'a', 'l', 'l', 'o', ' ', 'W', 'e', 'l', 't',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                                    };

    auto vector = Vector{text};

    REQUIRE(vector.size() == text.size());

    SECTION("check that symbol() call works") {
        for (size_t i{0}; i < text.size(); ++i) {
            INFO(i);
            CHECK(vector.symbol(i) == text.at(i));
        }
    }
    CHECK(Vector::Sigma >= 128);

    auto countRank = [&](size_t idx, uint8_t sym) {
        size_t acc{};
        for (size_t i{0}; i < idx; ++i) {
            acc = acc + (text[i] == sym);
        }
        return acc;
    };
    auto countPrefixRank = [&](size_t idx, uint8_t sym) {
        size_t acc{};
        for (size_t i{0}; i < idx; ++i) {
            acc = acc + (text[i] <= sym);
        }
        return acc;
    };


    SECTION("check all_ranks() is equal to prefix_rank() and rank()") {
        for (size_t idx{0}; idx <= vector.size(); ++idx) {
            auto [rank, prefix] = vector.all_ranks_and_prefix_ranks(idx);
            auto rank2 = vector.all_ranks(idx);
            for (size_t symb{1}; symb < Vector::Sigma; ++symb) {
                INFO(idx);
                INFO(symb);
                CHECK(countRank(idx, symb) == vector.rank(idx, symb));
                CHECK(countPrefixRank(idx, symb) == vector.prefix_rank(idx, symb));
                vector.rank(idx, symb);
                vector.prefix_rank(idx, symb);
                CHECK(rank[symb] == vector.rank(idx, symb));
                CHECK(rank2[symb] == vector.rank(idx, symb));
                CHECK(countRank(idx, symb) == rank[symb]);
                CHECK(countPrefixRank(idx, symb) == prefix[symb]);
                CHECK(prefix[symb] == vector.prefix_rank(idx, symb));
            }
        }
    }
}


static auto benchs_256 = Benchs{};
static auto benchSize_256 = BenchSize{};

TEMPLATE_TEST_CASE("benchmark vectors c'tor", "[RankVector][!benchmark][256][time][ctor][.]", ALLRANKVECTORS(256)) {
    using Vector = TestType;
    if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<256>>) {
        return;
    }

    auto& [bench_rank, bench_prefix_rank, bench_all_ranks, bench_all_prefix_ranks, bench_symbol, bench_ctor] = benchs_256;


    auto vector_name = getName<Vector>();
    INFO(vector_name);

    SECTION("benchmarking") {
        auto rng = ankerl::nanobench::Rng{};

        // generates string with values between 1-4
        auto text = std::vector<uint8_t>{};
        #ifdef NDEBUG
        for (size_t i{0}; i<1'000'000; ++i) {
        #else
        for (size_t i{0}; i<1'000; ++i) {
        #endif
            text.push_back(rng.bounded(256));
        }
        size_t minEpochIterations = 2'000'000;
        minEpochIterations = 1;

        bench_ctor.minEpochIterations(minEpochIterations).batch(text.size()).run(vector_name, [&]() {
            auto vec = Vector{text};
            ankerl::nanobench::doNotOptimizeAway(const_cast<Vector const&>(vec));
        });
    }
}

TEMPLATE_TEST_CASE("benchmark vectors symbol() and rank() operations", "[RankVector][!benchmark][256][time]", ALLRANKVECTORS(256)) {
    using Vector = TestType;
    if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<256>>) {
        return;
    }

    auto& [bench_rank, bench_prefix_rank, bench_all_ranks, bench_all_prefix_ranks, bench_symbol, bench_ctor] = benchs_256;


    auto vector_name = getName<Vector>();
    INFO(vector_name);

    SECTION("benchmarking") {
        auto rng = ankerl::nanobench::Rng{};

        // generates string with values between 1-4
        auto text = std::vector<uint8_t>{};
        #ifdef NDEBUG
        for (size_t i{0}; i<1'000'000; ++i) {
        #else
        for (size_t i{0}; i<1'000; ++i) {
        #endif
            text.push_back(rng.bounded(256));
        }
        auto vec = Vector{text};

        size_t minEpochIterations = 2'000'000;
        minEpochIterations = 1;

        bench_symbol.minEpochIterations(minEpochIterations).run(vector_name, [&]() {
            auto v = vec.symbol(rng.bounded(text.size()));
            ankerl::nanobench::doNotOptimizeAway(v);
        });

        bench_rank.minEpochIterations(minEpochIterations).run(vector_name, [&]() {
            auto v = vec.rank(rng.bounded(text.size()), rng.bounded(256));
            ankerl::nanobench::doNotOptimizeAway(v);
        });

        bench_prefix_rank.minEpochIterations(minEpochIterations).run(vector_name, [&]() {
            auto v = vec.prefix_rank(rng.bounded(text.size()), rng.bounded(256));
            ankerl::nanobench::doNotOptimizeAway(v);
        });

        bench_all_ranks.minEpochIterations(minEpochIterations).run(vector_name, [&]() {
            auto v = vec.all_ranks(rng.bounded(text.size()));
            ankerl::nanobench::doNotOptimizeAway(v);
        });

        bench_all_prefix_ranks.minEpochIterations(minEpochIterations).run(vector_name, [&]() {
            auto v = vec.all_ranks_and_prefix_ranks(rng.bounded(text.size()));
            ankerl::nanobench::doNotOptimizeAway(v);
        });
    }
}

TEMPLATE_TEST_CASE("benchmark vectors alphabet=256 size", "[RankVector][!benchmark][256][size]", ALLRANKVECTORS(256)) {
    using Vector = TestType;
    if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<256>>) {
        return;
    }

    auto vector_name = getName<Vector>();
    INFO(vector_name);

    SECTION("benchmarking") {
        auto rng = ankerl::nanobench::Rng{};

        // generates string with values between 1-4
        auto text = std::vector<uint8_t>{};
        #ifdef NDEBUG
        for (size_t i{0}; i<1'000'000; ++i) {
        #else
        for (size_t i{0}; i<1'000; ++i) {
        #endif
            text.push_back(rng.bounded(256));
        }
        auto vec = Vector{text};


        {
            auto ofs     = std::stringstream{};
            auto archive = cereal::BinaryOutputArchive{ofs};
            archive(vec);
            auto s = ofs.str().size();
            benchSize_256.addEntry({
                .name = vector_name,
                .size = s,
                .text_size = text.size(),
                .bits_per_char = (s*8)/double(text.size())
            });
        }
    }
}
#endif

#if 1
TEST_CASE("benchmark vectors c'tor operation, dna4 like", "[RankVector][!benchmark][5][time][ctor][.]") {
    auto const& text = generateText<1, 4>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("c'tor()")
             .relative(true)
             .batch(text.size());

        call_with_templates<
            ALLRANKVECTORS(5)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<5>>) {
                return;
            }

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            bench.run(vector_name, [&]() {
                auto vec = Vector{text};
                ankerl::nanobench::doNotOptimizeAway(const_cast<Vector const&>(vec));
            });
        });
    }
}

TEST_CASE("benchmark vectors symbol() operations, dna4 like", "[RankVector][!benchmark][5][time][symbol]") {
    auto const& text = generateText<1, 4>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("symbol()")
             .relative(true)
             .batch(text.size());

        call_with_templates<
            ALLRANKVECTORS(5)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<5>>) {
                return;
            }

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto rng = ankerl::nanobench::Rng{};

            auto vec = Vector{text};

            bench.run(vector_name, [&]() {
                auto v = vec.symbol(rng.bounded(text.size()));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        });
    }
}
#endif
TEST_CASE("benchmark vectors rank() operations, dna4 like", "[RankVector][!benchmark][5][time][rank]") {
    auto const& text = generateText<1, 4>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("rank()")
             .relative(true);

        call_with_templates<
            ALLRANKVECTORS(5)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<5>>) {
                return;
            }

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto rng = ankerl::nanobench::Rng{};

            auto vec = Vector{text};

            bench.run(vector_name, [&]() {
                auto v = vec.rank(rng.bounded(text.size()+1), rng.bounded(4)+1);
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        });
    }
}
#if 1
TEST_CASE("benchmark vectors prefix_rank() operations, dna4 like", "[RankVector][!benchmark][5][time][prefix_rank]") {
    auto const& text = generateText<1, 4>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("prefix_rank()")
             .relative(true);

        call_with_templates<
            ALLRANKVECTORS(5)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<5>>) {
                return;
            }

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto rng = ankerl::nanobench::Rng{};

            auto vec = Vector{text};

            bench.run(vector_name, [&]() {
                auto v = vec.prefix_rank(rng.bounded(text.size()+1), rng.bounded(4)+1);
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        });
    }
}

TEST_CASE("benchmark vectors all_ranks() operations, dna4 like", "[RankVector][!benchmark][5][time][all_ranks]") {
    auto const& text = generateText<1, 4>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("all_ranks()")
             .relative(true);

        call_with_templates<
            ALLRANKVECTORS(5)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<5>>) {
                return;
            }

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto rng = ankerl::nanobench::Rng{};

            auto vec = Vector{text};

            bench.run(vector_name, [&]() {
                auto v = vec.all_ranks(rng.bounded(text.size()+1));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        });
    }
}
TEST_CASE("benchmark vectors all_ranks_and_prefix_ranks() operations, dna4 like", "[RankVector][!benchmark][5][time][all_ranks_and_prefix_ranks]") {
    auto const& text = generateText<1, 4>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("all_ranks_and_prefix_ranks()")
             .relative(true);

        call_with_templates<
            ALLRANKVECTORS(5)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<5>>) {
                return;
            }

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto rng = ankerl::nanobench::Rng{};

            auto vec = Vector{text};

            bench.run(vector_name, [&]() {
                auto v = vec.all_ranks_and_prefix_ranks(rng.bounded(text.size()+1));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        });
    }
}

TEST_CASE("benchmark vectors in size, dna4 like", "[RankVector][!benchmark][5][size]") {
    auto const& text = generateText<1, 4>();

    SECTION("benchmarking") {
        BenchSize benchSize;
        benchSize.entries[0][2] = "bits/char";

        call_with_templates<
            ALLRANKVECTORS(5)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<5>>) {
                return;
            }

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto rng = ankerl::nanobench::Rng{};

            auto vec = Vector{text};

            {
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
#endif
/*
TEMPLATE_TEST_CASE("benchmark vectors alphabet=5 (dna4 like) in size", "[RankVector][!benchmark][5][size]", ALLRANKVECTORS(5)) {
    using Vector = TestType;
    if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<5>>) {
        return;
    }

    auto vector_name = getName<Vector>();
    INFO(vector_name);

    SECTION("benchmarking") {
        auto rng = ankerl::nanobench::Rng{};

        // generates string with values between 1-4
        auto text = std::vector<uint8_t>{};
        #ifdef NDEBUG
        for (size_t i{0}; i<1'000'000; ++i) {
        #else
        for (size_t i{0}; i<1'000; ++i) {
        #endif
            text.push_back(rng.bounded(4)+1);
        }
        auto vec = Vector{text};

        {
            auto ofs     = std::stringstream{};
            auto archive = cereal::BinaryOutputArchive{ofs};
            archive(vec);
            auto s = ofs.str().size();
            benchSize_5.addEntry({
                .name = vector_name,
                .size = s,
                .text_size = text.size(),
                .bits_per_char = (s*8)/double(text.size())
            });
        }
    }
}*/
/*
    using Vector = TestType;
    if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<5>>) {
        return;
    }
    auto& [bench_rank, bench_prefix_rank, bench_all_ranks, bench_all_prefix_ranks, bench_symbol, bench_ctor] = benchs_5;


    auto vector_name = getName<Vector>();
    INFO(vector_name);

    SECTION("benchmarking") {
        auto rng = ankerl::nanobench::Rng{};

        // generates string with values between 1-4
        auto text = std::vector<uint8_t>{};
        #ifdef NDEBUG
        for (size_t i{0}; i<1'000'000; ++i) {
        #else
        for (size_t i{0}; i<1'000; ++i) {
        #endif
            text.push_back(rng.bounded(4)+1);
        }
        auto vec = Vector{text};

        size_t minEpochIterations = 2'000'000;
        minEpochIterations = 1;


        bench_symbol.minEpochIterations(minEpochIterations).run(vector_name, [&]() {
            auto v = vec.symbol(rng.bounded(text.size()));
            ankerl::nanobench::doNotOptimizeAway(v);
        });

        bench_rank.minEpochIterations(minEpochIterations).run(vector_name, [&]() {
            auto v = vec.rank(rng.bounded(text.size()), rng.bounded(4)+1);
            ankerl::nanobench::doNotOptimizeAway(v);
        });

        bench_prefix_rank.minEpochIterations(minEpochIterations).run(vector_name, [&]() {
            auto v = vec.prefix_rank(rng.bounded(text.size()), rng.bounded(4)+1);
            ankerl::nanobench::doNotOptimizeAway(v);
        });

        bench_all_ranks.minEpochIterations(minEpochIterations).run(vector_name, [&]() {
            auto v = vec.all_ranks(rng.bounded(text.size()));
            ankerl::nanobench::doNotOptimizeAway(v);
        });

        bench_all_prefix_ranks.minEpochIterations(minEpochIterations).run(vector_name, [&]() {
            auto v = vec.all_ranks_and_prefix_ranks(rng.bounded(text.size()));
            ankerl::nanobench::doNotOptimizeAway(v);
        });
    }
}
*/

#if 1
static auto benchs_6 = Benchs{};
static auto benchs_6_text = std::vector<uint8_t>{};
static auto benchSize_bwt = BenchSize{};
TEMPLATE_TEST_CASE("benchmark vectors c'tor operation, on human dna5 data", "[RankVector][bwt][!benchmark][6][time][ctor][.]", ALLRANKVECTORS(6)) {
    using Vector = TestType;
    if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<6>>) {
        return;
    }
    auto& [bench_rank, bench_prefix_rank, bench_all_ranks, bench_all_prefix_ranks, bench_symbol, bench_ctor] = benchs_6;


    auto vector_name = getName<Vector>();
    INFO(vector_name);

    SECTION("benchmarking") {
        auto rng = ankerl::nanobench::Rng{};

        // generates string with values between 1-4
        auto& text = benchs_6_text;
        if (text.empty()) {
            auto fs = std::ifstream{"bwt.bwt"};
            std::string line;
            std::getline(fs, line);
            fs.close();
            text.resize(line.size());
            std::memcpy(text.data(), line.data(), line.size());
            for (auto& c : text) {
                switch(c) {
                case '$': c = 0; break;
                case 'A': c = 1; break;
                case 'C': c = 2; break;
                case 'G': c = 3; break;
                case 'T': c = 4; break;
                case 'N': c = 5; break;
                default:
                    throw std::runtime_error{"unknown char " + std::to_string((int)c)};

                }
            }
        }

        auto vec = Vector{text};


        size_t minEpochIterations = 2'000'000;
        minEpochIterations = 1;
        bench_ctor.minEpochIterations(minEpochIterations).batch(text.size()).run(vector_name, [&]() {
            auto vec = Vector{text};
            ankerl::nanobench::doNotOptimizeAway(const_cast<Vector const&>(vec));
        });
    }
}

TEMPLATE_TEST_CASE("benchmark vectors symbol() and rank() operations, on human dna5 data", "[RankVector][bwt][!benchmark][6][time][.]", ALLRANKVECTORS(6)) {
    using Vector = TestType;
    if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<6>>) {
        return;
    }
    auto& [bench_rank, bench_prefix_rank, bench_all_ranks, bench_all_prefix_ranks, bench_symbol, bench_ctor] = benchs_6;


    auto vector_name = getName<Vector>();
    INFO(vector_name);

    SECTION("benchmarking") {
        auto rng = ankerl::nanobench::Rng{};

        // generates string with values between 1-4
        auto& text = benchs_6_text;
        if (text.empty()) {
            auto fs = std::ifstream{"bwt.bwt"};
            std::string line;
            std::getline(fs, line);
            fs.close();
            text.resize(line.size());
            std::memcpy(text.data(), line.data(), line.size());
            for (auto& c : text) {
                switch(c) {
                case '$': c = 0; break;
                case 'A': c = 1; break;
                case 'C': c = 2; break;
                case 'G': c = 3; break;
                case 'T': c = 4; break;
                case 'N': c = 5; break;
                default:
                    throw std::runtime_error{"unknown char " + std::to_string((int)c)};

                }
            }
        }

        auto vec = Vector{text};
        REQUIRE(text.size() > 0);


        size_t minEpochIterations = 2'000'000;
        minEpochIterations = 1;
        bench_symbol.minEpochIterations(minEpochIterations).run(vector_name, [&]() {
            auto v = vec.symbol(rng.bounded(text.size()));
            ankerl::nanobench::doNotOptimizeAway(v);
        });

        bench_rank.minEpochIterations(minEpochIterations).run(vector_name, [&]() {
            auto v = vec.rank(rng.bounded(text.size()), rng.bounded(5)+1);
            ankerl::nanobench::doNotOptimizeAway(v);
        });

        bench_prefix_rank.minEpochIterations(minEpochIterations).run(vector_name, [&]() {
            auto v = vec.prefix_rank(rng.bounded(text.size()), rng.bounded(5)+1);
            ankerl::nanobench::doNotOptimizeAway(v);
        });

        bench_all_ranks.minEpochIterations(minEpochIterations).run(vector_name, [&]() {
            auto v = vec.all_ranks(rng.bounded(text.size()));
            ankerl::nanobench::doNotOptimizeAway(v);
        });

        bench_all_prefix_ranks.minEpochIterations(minEpochIterations).run(vector_name, [&]() {
            auto v = vec.all_ranks_and_prefix_ranks(rng.bounded(text.size()));
            ankerl::nanobench::doNotOptimizeAway(v);
        });
    }
}

TEMPLATE_TEST_CASE("benchmark vectors size, on human dna5 data", "[RankVector][bwt][!benchmark][6][size][.]", ALLRANKVECTORS(6)) {
    using Vector = TestType;
    if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<6>>) {
        return;
    }

    auto vector_name = getName<Vector>();
    INFO(vector_name);

    SECTION("benchmarking") {
        auto rng = ankerl::nanobench::Rng{};

        // generates string with values between 1-4
        auto& text = benchs_6_text;
        if (text.empty()) {
            auto fs = std::ifstream{"bwt.bwt"};
            std::string line;
            std::getline(fs, line);
            fs.close();
            text.resize(line.size());
            std::memcpy(text.data(), line.data(), line.size());
            for (auto& c : text) {
                switch(c) {
                case '$': c = 0; break;
                case 'A': c = 1; break;
                case 'C': c = 2; break;
                case 'G': c = 3; break;
                case 'T': c = 4; break;
                case 'N': c = 5; break;
                default:
                    throw std::runtime_error{"unknown char " + std::to_string((int)c)};

                }
            }
        }

        auto vec = Vector{text};

        auto ofs     = std::stringstream{};
        auto archive = cereal::BinaryOutputArchive{ofs};
        archive(vec);
        auto s = ofs.str().size();
        benchSize_bwt.addEntry({
            .name = vector_name,
            .size = s,
            .text_size = text.size(),
            .bits_per_char = (s*8)/double(text.size())
        });
    }
}

TEST_CASE("benchmark vectors c'tor operation, 255 like", "[RankVector][!benchmark][255][time][ctor][.]") {
    auto const& text = generateText<0, 255>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("c'tor()")
             .relative(true)
             .batch(text.size());

        call_with_templates<
            ALLRANKVECTORS(255)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<255>>) {
                return;
            }

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            bench.run(vector_name, [&]() {
                auto vec = Vector{text};
                ankerl::nanobench::doNotOptimizeAway(const_cast<Vector const&>(vec));
            });
        });
    }
}

TEST_CASE("benchmark vectors symbol() operations, 255 like", "[RankVector][!benchmark][255][time][symbol]") {
    auto const& text = generateText<0, 255>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("symbol()")
             .relative(true)
             .batch(text.size());

        call_with_templates<
            ALLRANKVECTORS(255)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<255>>) {
                return;
            }

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto rng = ankerl::nanobench::Rng{};

            auto vec = Vector{text};

            bench.run(vector_name, [&]() {
                auto v = vec.symbol(rng.bounded(text.size()));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        });
    }
}
TEST_CASE("benchmark vectors rank() operations, 255 like", "[RankVector][!benchmark][255][time][rank]") {
    auto const& text = generateText();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("rank()")
             .relative(true);

        call_with_templates<
            ALLRANKVECTORS(255)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<255>>) {
                return;
            }

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto rng = ankerl::nanobench::Rng{};

            auto vec = Vector{text};

            bench.run(vector_name, [&]() {
                auto v = vec.rank(rng.bounded(text.size()+1), rng.bounded(4)+1);
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        });
    }
}

TEST_CASE("benchmark vectors prefix_rank() operations, 255 like", "[RankVector][!benchmark][255][time][prefix_rank]") {
    auto const& text = generateText<0, 255>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("prefix_rank()")
             .relative(true);

        call_with_templates<
            ALLRANKVECTORS(255)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<255>>) {
                return;
            }

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto rng = ankerl::nanobench::Rng{};

            auto vec = Vector{text};

            bench.run(vector_name, [&]() {
                auto v = vec.prefix_rank(rng.bounded(text.size()+1), rng.bounded(4)+1);
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        });
    }
}

TEST_CASE("benchmark vectors all_ranks() operations, 255 like", "[RankVector][!benchmark][255][time][all_ranks]") {
    auto const& text = generateText<0, 255>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("all_ranks()")
             .relative(true);

        call_with_templates<
            ALLRANKVECTORS(255)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<255>>) {
                return;
            }

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto rng = ankerl::nanobench::Rng{};

            auto vec = Vector{text};

            bench.run(vector_name, [&]() {
                auto v = vec.all_ranks(rng.bounded(text.size()+1));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        });
    }
}
TEST_CASE("benchmark vectors all_ranks_and_prefix_ranks() operations, 255 like", "[RankVector][!benchmark][255][time][all_ranks_and_prefix_ranks]") {
    auto const& text = generateText<0, 255>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("all_ranks_and_prefix_ranks()")
             .relative(true);

        call_with_templates<
            ALLRANKVECTORS(255)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<255>>) {
                return;
            }

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto rng = ankerl::nanobench::Rng{};

            auto vec = Vector{text};

            bench.run(vector_name, [&]() {
                auto v = vec.all_ranks_and_prefix_ranks(rng.bounded(text.size()+1));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        });
    }
}

TEST_CASE("benchmark vectors in size, 255 like", "[RankVector][!benchmark][255][size]") {
    auto const& text = generateText<0, 255>();

    SECTION("benchmarking") {
        BenchSize benchSize;
        benchSize.entries[0][2] = "bits/char";

        call_with_templates<
            ALLRANKVECTORS(255)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<255>>) {
                return;
            }

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto rng = ankerl::nanobench::Rng{};

            auto vec = Vector{text};

            {
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

#endif
