// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "utils.h"

TEST_CASE("benchmark strings c'tor operation - 4 alphabet", "[string][!benchmark][4][time][ctor]") {
    auto const& text = generateText<0, 4>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("c'tor()")
             .relative(true)
             .batch(text.size());

        call_with_templates<
            ALLRANKVECTORS(4)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::string::Naive<4>>) {
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

TEST_CASE("benchmark vectors symbol() operations - 4 alphabet", "[string][!benchmark][4][time][symbol]") {
    auto const& text = generateText<0, 4>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("symbol()")
             .relative(true)
             .batch(text.size());

        call_with_templates<
            ALLRANKVECTORS(4)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::string::Naive<4>>) {
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

TEST_CASE("benchmark vectors rank() operations - 4 alphabet", "[string][!benchmark][4][time][rank]") {
    auto const& text = generateText<0, 4>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("rank()")
             .relative(true);

        call_with_templates<
            ALLRANKVECTORS(4)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::string::Naive<4>>) {
                return;
            }

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto rng = ankerl::nanobench::Rng{};

            auto vec = Vector{text};

            bench.run(vector_name, [&]() {
                auto v = vec.rank(rng.bounded(text.size()+1), rng.bounded(4));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        });
    }
}

TEST_CASE("benchmark vectors prefix_rank() operations - 4 alphabet", "[string][!benchmark][4][time][prefix_rank]") {
    auto const& text = generateText<0, 4>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("prefix_rank()")
             .relative(true);

        call_with_templates<
            ALLRANKVECTORS(4)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::string::Naive<4>>) {
                return;
            }

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto rng = ankerl::nanobench::Rng{};

            auto vec = Vector{text};

            bench.run(vector_name, [&]() {
                auto v = vec.prefix_rank(rng.bounded(text.size()+1), rng.bounded(4));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        });
    }
}

TEST_CASE("benchmark vectors all_ranks() operations - 4 alphabet", "[string][!benchmark][4][time][all_ranks]") {
    auto const& text = generateText<0, 4>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("all_ranks()")
             .relative(true);

        call_with_templates<
            ALLRANKVECTORS(4)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::string::Naive<4>>) {
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

TEST_CASE("benchmark vectors all_ranks_and_prefix_ranks() operations - 4 alphabet", "[string][!benchmark][4][time][all_ranks_and_prefix_ranks]") {
    auto const& text = generateText<0, 4>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("all_ranks_and_prefix_ranks()")
             .relative(true);

        call_with_templates<
            ALLRANKVECTORS(4)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::string::Naive<4>>) {
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

TEST_CASE("benchmark vectors in size - alphabet 4", "[string][!benchmark][4][size]") {
    auto const& text = generateText<0, 4>();

    SECTION("benchmarking") {
        BenchSize benchSize;
        benchSize.entries[0][2] = "bits/char";
        benchSize.entries[0][3] = "alphabet 4";

        call_with_templates<
            ALLRANKVECTORS(4)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::string::Naive<4>>) {
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
