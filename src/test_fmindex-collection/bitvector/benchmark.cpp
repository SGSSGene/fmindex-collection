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

#ifdef FMC_USE_PASTA
#include "Pasta_FlatRank.h"
#include "Pasta_WideRank.h"
#endif

#ifdef FMC_USE_SDSL
#include "sdsl_v.h"
#include "sdsl_v5.h"
#endif

#ifdef FMC_USE_SUX
#include "sux_Rank9.h"
#endif

#ifdef FMC_USE_RANKSELECT
#include "RankSelect_Rank.h"
#endif

namespace {
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
#if defined(FMC_USE_PASTA) && defined(FMC_USE_SDSL) && defined(FMC_USE_SUX) && defined(FMC_USE_RANKSELECT)
    #define ALLTYPES \
             ALLBITVECTORS, \
             FlatRank, \
             WideRank, \
             SDSL_V, \
             SDSL_V5, \
             Rank9, \
             RankSelect<0>

#else
    #define ALLTYPES ALLBITVECTORS
#endif


TEST_CASE("benchmark bit vectors ctor run times", "[bitvector][!benchmark][time][ctor]") {
    auto bench_ctor = ankerl::nanobench::Bench{};
    bench_ctor.title("c'tor()")
              .relative(true);

    auto& text = generateText();

    SECTION("benchmarking") {
        call_with_templates<ALLTYPES>([&]<typename Vector>() {

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            bench_ctor.batch(text.size()).run(vector_name, [&]() {
                auto vec = Vector{text};
                ankerl::nanobench::doNotOptimizeAway(vec.rank(0));
            });
        });
    }
}

TEST_CASE("benchmark bit vectors rank and symbol run times", "[bitvector][!benchmark][time][symbol]") {

    auto& text = generateText();

    SECTION("benchmarking - symbol") {
        auto bench_symbol = ankerl::nanobench::Bench{};
        bench_symbol.title("symbol()")
                    .relative(true);

        bench_symbol.epochs(10);
        bench_symbol.minEpochTime(std::chrono::milliseconds{10});
        call_with_templates<ALLTYPES>([&]<typename Vector>() {

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

TEST_CASE("benchmark bit vectors rank and symbol run times", "[bitvector][!benchmark][time][rank]") {

    auto& text = generateText();

    SECTION("benchmarking - rank") {
        auto bench_rank = ankerl::nanobench::Bench{};
        bench_rank.title("rank()")
                  .relative(true);

        bench_rank.epochs(20);
        bench_rank.minEpochTime(std::chrono::milliseconds{1});
        bench_rank.minEpochIterations(1'000'000);

        call_with_templates<ALLTYPES>([&]<typename Vector>() {

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

TEST_CASE("benchmark bit vectors memory consumption", "[bitvector][!benchmark][size]") {
    BenchSize benchSize;
    benchSize.baseSize = 1.;

    SECTION("benchmarking") {
        call_with_templates<ALLTYPES>([&]<typename Vector>() {

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
