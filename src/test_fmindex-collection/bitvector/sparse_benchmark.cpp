// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include <catch2/catch_all.hpp>
#include <cstdlib>
#include <fmindex-collection/utils.h>
#include <fstream>
#include <nanobench.h>

#include "../string/utils.h"
#include "../BenchSize.h"
#include "allBitVectors.h"

using AllTypes = std::variant<
    ALLSPARSEBITVECTORS,
#if defined(FMC_USE_RANKSELECT)
    RankSelect<0>,
    RankSelect<1>,
    RankSelect<2>,
    RankSelect<3>,
    RankSelect<4>,
    RankSelect<5>,
    RankSelect<6>,
    RankSelect<7>,
#endif
    std::monostate
>;

namespace {
auto generateText(double density) -> std::vector<bool> {
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
        text.push_back(rng.bounded(100'000) < density*100'000.);
    }
    return text;
}

auto generateText(double density, size_t blockSize) -> std::vector<bool> {
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
    for (size_t i{0}; i+blockSize-1 < size; i += blockSize) {
        auto v = (rng.bounded(100'000) < density*100'000);

        for (size_t j{0}; j < blockSize; ++j) {
            text.push_back(v);
        }
    }
    while (text.size() < size) {
        text.push_back(rng.bounded(100'000) < density*100'000);
    }
    return text;
}

}


TEST_CASE("benchmark bit vectors ctor run times", "[sparse-bitvector][!benchmark][time][ctor]") {
    auto bench_ctor = ankerl::nanobench::Bench{};
    bench_ctor.title("c'tor()")
              .relative(true);

    auto text = generateText(0.5);

    SECTION("benchmarking") {
        call_with_templates<AllTypes>([&]<typename Vector>() {

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            bench_ctor.batch(text.size()).run(vector_name, [&]() {
                auto vec = Vector{text};
                ankerl::nanobench::doNotOptimizeAway(vec.rank(0));
            });
        });
    }
}

TEST_CASE("benchmark bit vectors rank and symbol run times", "[sparse-bitvector][!benchmark][time][symbol]") {

    auto text = generateText(0.5);

    SECTION("benchmarking - symbol") {
        auto bench_symbol = ankerl::nanobench::Bench{};
        bench_symbol.title("symbol()")
                    .relative(true);

        bench_symbol.epochs(10);
        bench_symbol.minEpochTime(std::chrono::milliseconds{10});
        call_with_templates<AllTypes>([&]<typename Vector>() {

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
}

TEST_CASE("benchmark bit vectors rank and symbol run times", "[sparse-bitvector][!benchmark][time][rank]") {

    SECTION("benchmarking - rank") {
        for (auto density : {0.05, 0.10, 0.15, 0.20, 0.25, 0.3, 0.35, 0.40, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95}) {
            auto text = generateText(density);

            auto bench_rank = ankerl::nanobench::Bench{};
            bench_rank.title("rank() density:" + std::to_string(density))
                      .relative(true);

            bench_rank.epochs(20);
            bench_rank.minEpochTime(std::chrono::milliseconds{1});
            bench_rank.minEpochIterations(1'000'000);

            call_with_templates<AllTypes>([&]<typename Vector>() {

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
}

TEST_CASE("benchmark bit vectors memory consumption", "[sparse-bitvector][!benchmark][size][density]") {
    SECTION("benchmarking") {
        for (auto density : {0.05, 0.10, 0.15, 0.20, 0.25, 0.3, 0.35, 0.40, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95}) {
            BenchSize benchSize;
            benchSize.baseSize = 1.;

            auto text = generateText(density);


            call_with_templates<AllTypes>([&]<typename Vector>() {

                auto vector_name = getName<Vector>() + " density " + std::to_string(density);
                INFO(vector_name);

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
}

TEST_CASE("benchmark bit vectors memory consumption - with blocks", "[sparse-bitvector][!benchmark][size][block]") {
    SECTION("benchmarking") {
        for (auto density : {0.05, 0.10, 0.15, 0.20, 0.25, 0.3, 0.35, 0.40, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95}) {
            for (size_t blockiness{2}; blockiness < 3; ++blockiness) {
//            for (size_t blockiness{0}; blockiness < 6; ++blockiness) {
                BenchSize benchSize;
                benchSize.baseSize = 1.;

                auto blockSize = 1<<blockiness;

                auto text = generateText(density, blockSize);

                call_with_templates<AllTypes>([&]<typename Vector>() {

                    auto vector_name = getName<Vector>() + " density " + std::to_string(density) + ", blocksize " + std::to_string(blockSize);
                    INFO(vector_name);

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
    }
}
