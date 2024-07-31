// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include <catch2/catch_all.hpp>
#include <fmindex-collection/utils.h>
#include <nanobench.h>
#include <reflect>

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
template <typename T>
static auto getName() {
    return std::string{reflect::type_name<T>()};
}
#endif

TEMPLATE_TEST_CASE("check if rank on the symbol vectors is working", "[RankVector]", ALLRANKVECTORS(256)) {
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
    }

    SECTION("test complete vec 'H' for prefix_rank()") {
        CHECK(vec.prefix_rank( 0, ' ') == 0);
        CHECK(vec.prefix_rank( 1, ' ') == 0);
        CHECK(vec.prefix_rank( 2, ' ') == 0);
        CHECK(vec.prefix_rank( 3, ' ') == 0);
        CHECK(vec.prefix_rank( 4, ' ') == 0);
        CHECK(vec.prefix_rank( 5, ' ') == 0);
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

namespace {
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

static auto benchs_256 = Benchs{};
TEMPLATE_TEST_CASE("benchmark vectors c'tor,symbol() and rank() operations", "[RankVector][!benchmark]", ALLRANKVECTORS(256)) {
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


        bench_ctor.minEpochIterations(10).batch(text.size()).run(vector_name, [&]() {
            auto vec = Vector{text};
            ankerl::nanobench::doNotOptimizeAway(const_cast<Vector const&>(vec));
        });

        size_t minEpochIterations = 2'000'000;
        minEpochIterations = 2'000;
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

        {
            auto ofs     = std::stringstream{};
            auto archive = cereal::BinaryOutputArchive{ofs};
            archive(vec);
            auto s = ofs.str().size()*8;
            std::cout << vector_name << " - file size: " << s << "bytes, " << s/double(text.size()) << "bits/bits\n";
        }
    }
}

static auto benchs_5 = Benchs{};
TEMPLATE_TEST_CASE("benchmark vectors c'tor,symbol() and rank() operations, dna4 like", "[RankVector][!benchmark]", ALLRANKVECTORS(5)) {
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


        bench_ctor.minEpochIterations(10).batch(text.size()).run(vector_name, [&]() {
            auto vec = Vector{text};
            ankerl::nanobench::doNotOptimizeAway(const_cast<Vector const&>(vec));
        });

        bench_symbol.minEpochIterations(2'000'000).run(vector_name, [&]() {
            auto v = vec.symbol(rng.bounded(text.size()));
            ankerl::nanobench::doNotOptimizeAway(v);
        });

        bench_rank.minEpochIterations(2'000'000).run(vector_name, [&]() {
            auto v = vec.rank(rng.bounded(text.size()), rng.bounded(4)+1);
            ankerl::nanobench::doNotOptimizeAway(v);
        });

        bench_prefix_rank.minEpochIterations(2'000'000).run(vector_name, [&]() {
            auto v = vec.prefix_rank(rng.bounded(text.size()), rng.bounded(4)+1);
            ankerl::nanobench::doNotOptimizeAway(v);
        });

        bench_all_ranks.minEpochIterations(2'000'000).run(vector_name, [&]() {
            auto v = vec.all_ranks(rng.bounded(text.size()));
            ankerl::nanobench::doNotOptimizeAway(v);
        });

        bench_all_prefix_ranks.minEpochIterations(2'000'000).run(vector_name, [&]() {
            auto v = vec.all_ranks_and_prefix_ranks(rng.bounded(text.size()));
            ankerl::nanobench::doNotOptimizeAway(v);
        });

        {
            auto ofs     = std::stringstream{};
            auto archive = cereal::BinaryOutputArchive{ofs};
            archive(vec);
            auto s = ofs.str().size()*8;
            std::cout << vector_name << " - file size: " << s << "bytes, " << s/double(text.size()) << "bits/bits\n";
        }
    }
}
