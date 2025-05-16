// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "utils.h"

TEST_CASE("benchmark strings c'tor operation - 255 alphabet", "[string][!benchmark][255][time][ctor]") {
    auto const& text = generateText<0, 255>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("c'tor()")
             .relative(true)
             .batch(text.size());

        call_with_templates<
            STRINGSWITHRANK(255)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::string::Naive<255>>) {
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

TEST_CASE("benchmark vectors symbol() operations - 255 alphabet", "[string][!benchmark][255][time][symbol]") {
    auto const& text = generateText<0, 255>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("symbol()")
             .relative(true)
             .batch(text.size());

        call_with_templates<
            STRINGSWITHRANK(255)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::string::Naive<255>>) {
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

TEST_CASE("benchmark vectors rank() operations - 255 alphabet", "[string][!benchmark][255][time][rank]") {
    auto const& text = generateText<0, 255>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("rank()")
             .relative(true);

        call_with_templates<
            STRINGSWITHRANK(255)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::string::Naive<255>>) {
                return;
            }

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto rng = ankerl::nanobench::Rng{};

            auto vec = Vector{text};

            bench.run(vector_name, [&]() {
                auto v = vec.rank(rng.bounded(text.size()+1), rng.bounded(255));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        });
    }
}

TEST_CASE("benchmark vectors prefix_rank() operations - 255 alphabet", "[string][!benchmark][255][time][prefix_rank]") {
    auto const& text = generateText<0, 255>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("prefix_rank()")
             .relative(true);

        call_with_templates<
            STRINGSWITHRANK(255)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::string::Naive<255>>) {
                return;
            }

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto rng = ankerl::nanobench::Rng{};

            auto vec = Vector{text};

            bench.run(vector_name, [&]() {
                auto v = vec.prefix_rank(rng.bounded(text.size()+1), rng.bounded(255));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        });
    }
}

TEST_CASE("benchmark vectors all_ranks() operations - 255 alphabet", "[string][!benchmark][255][time][all_ranks]") {
    auto const& text = generateText<0, 255>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("all_ranks()")
             .relative(true);

        call_with_templates<
            STRINGSWITHRANK(255)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::string::Naive<255>>) {
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

TEST_CASE("benchmark vectors all_ranks_and_prefix_ranks() operations - 255 alphabet", "[string][!benchmark][255][time][all_ranks_and_prefix_ranks]") {
    auto const& text = generateText<0, 255>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("all_ranks_and_prefix_ranks()")
             .relative(true);

        call_with_templates<
            STRINGSWITHRANK(255)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::string::Naive<255>>) {
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

TEST_CASE("benchmark vectors in size - alphabet 255", "[string][!benchmark][255][size]") {
    auto const& text = generateText<0, 255>();

    SECTION("benchmarking") {
        BenchSize benchSize;
        benchSize.baseSize = 8.;
        benchSize.entries[0][2] = "bits/char";
        benchSize.entries[0][4] = "alphabet 255";

        call_with_templates<
            STRINGSWITHRANK(255)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::string::Naive<255>>) {
                return;
            }

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto rng = ankerl::nanobench::Rng{};

            auto vec = Vector{text};
            auto size = [&]() {
                if constexpr (requires() { vec.space_usage(); }) {
                    return vec.space_usage();
                } else {
                    auto ofs     = std::stringstream{};
                    auto archive = cereal::BinaryOutputArchive{ofs};
                    archive(vec);
                    return ofs.str().size();
                }
            }();
            benchSize.addEntry({
                .name = vector_name,
                .size = size,
                .text_size = text.size(),
                .bits_per_char = (size*8)/double(text.size())
            });
        });
    }
}
