// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "utils.h"

namespace {
    constexpr static size_t Sigma = 16384;
    #define SIGMA_STR "16384"
    #define ALLSTRINGS ALLLARGESTRINGSWITHRANK
}

TEST_CASE("benchmark strings c'tor operation - " SIGMA_STR " alphabet", "[string][!benchmark][" SIGMA_STR "][time][ctor]") {
    auto const& text = generateLargeText<0, Sigma>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("c'tor()")
             .relative(true)
             .batch(text.size());

        call_with_templates<ALLSTRINGS>([&]<template <size_t> class _String>() {
            using String = _String<Sigma>;
            auto name = getName<String>();
            INFO(name);

            bench.run(name, [&]() {
                auto str = String{text};
                ankerl::nanobench::doNotOptimizeAway(const_cast<String const&>(str));
            });
        });
    }
}

TEST_CASE("benchmark vectors symbol() operations - " SIGMA_STR " alphabet", "[string][!benchmark][" SIGMA_STR "][time][symbol]") {
    auto const& text = generateLargeText<0, Sigma>();
    auto rng = ankerl::nanobench::Rng{};

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("symbol()")
             .relative(true)
             .batch(text.size());

        call_with_templates<ALLSTRINGS>([&]<template <size_t> class _String>() {
            using String = _String<Sigma>;
            auto name = getName<String>();
            INFO(name);

            auto str = String{text};

            bench.run(name, [&]() {
                auto v = str.symbol(rng.bounded(text.size()));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        });
    }
}

TEST_CASE("benchmark vectors rank() operations - " SIGMA_STR " alphabet", "[string][!benchmark][" SIGMA_STR "][time][rank]") {
    auto const& text = generateLargeText<0, Sigma>();
    auto rng = ankerl::nanobench::Rng{};

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("rank()")
             .relative(true);

        call_with_templates<ALLSTRINGS>([&]<template <size_t> class _String>() {
            using String = _String<Sigma>;
            auto name = getName<String>();
            INFO(name);

            auto str = String{text};

            bench.run(name, [&]() {
                auto v = str.rank(rng.bounded(text.size()+1), rng.bounded(Sigma));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        });
    }
}

TEST_CASE("benchmark vectors prefix_rank() operations - " SIGMA_STR " alphabet", "[string][!benchmark][" SIGMA_STR "][time][prefix_rank]") {
    auto const& text = generateLargeText<0, Sigma>();
    auto rng = ankerl::nanobench::Rng{};

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("prefix_rank()")
             .relative(true);

        call_with_templates<ALLSTRINGS>([&]<template <size_t> class _String>() {
            using String = _String<Sigma>;
            auto name = getName<String>();
            INFO(name);

            auto str = String{text};

            bench.run(name, [&]() {
                auto v = str.prefix_rank(rng.bounded(text.size()+1), rng.bounded(Sigma));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        });
    }
}

TEST_CASE("benchmark vectors all_ranks() operations - " SIGMA_STR " alphabet", "[string][!benchmark][" SIGMA_STR "][time][all_ranks]") {
    auto const& text = generateLargeText<0, Sigma>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("all_ranks()")
             .relative(true);

        call_with_templates<ALLSTRINGS>([&]<template <size_t> class _String>() {
            using String = _String<Sigma>;
            auto name = getName<String>();
            INFO(name);

            auto rng = ankerl::nanobench::Rng{};

            auto str = String{text};

            bench.run(name, [&]() {
                auto v = str.all_ranks(rng.bounded(text.size()+1));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        });
    }
}

TEST_CASE("benchmark vectors all_ranks_and_prefix_ranks() operations - Sigma alphabet", "[string][!benchmark][" SIGMA_STR "][time][all_ranks_and_prefix_ranks]") {
    auto const& text = generateLargeText<0, Sigma>();
    auto rng = ankerl::nanobench::Rng{};

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("all_ranks_and_prefix_ranks()")
             .relative(true);

        call_with_templates<ALLSTRINGS>([&]<template <size_t> class _String>() {
            using String = _String<Sigma>;
            auto name = getName<String>();
            INFO(name);

            auto str = String{text};

            bench.run(name, [&]() {
                auto v = str.all_ranks_and_prefix_ranks(rng.bounded(text.size()+1));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        });
    }
}

TEST_CASE("benchmark vectors in size - alphabet " SIGMA_STR, "[string][!benchmark][" SIGMA_STR "][size]") {
    auto const& text = generateLargeText<0, Sigma>();
    auto rng = ankerl::nanobench::Rng{};

    SECTION("benchmarking") {
        BenchSize benchSize;
        benchSize.baseSize = 14.;
        benchSize.entries[0][2] = "bits/char";
        benchSize.entries[0][4] = "alphabet " SIGMA_STR;

        call_with_templates<ALLSTRINGS>([&]<template <size_t> class _String>() {
            using String = _String<Sigma>;
            auto name = getName<String>();
            INFO(name);

            auto str = String{text};
            auto size = [&]() {
                auto ofs     = std::stringstream{};
                auto archive = cereal::BinaryOutputArchive{ofs};
                archive(str);
                return ofs.str().size();
            }();
            benchSize.addEntry({
                .name = name,
                .size = size,
                .text_size = text.size(),
                .bits_per_char = (size*8)/double(text.size())
            });
        });
    }
}
