#include "utils.h"

TEST_CASE("benchmark strings c'tor operation - 128 alphabet", "[string_with_rank][!benchmark][128][time][ctor][.]") {
    auto const& text = generateText<0, 128>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("c'tor()")
             .relative(true)
             .batch(text.size());

        call_with_templates<
            ALLRANKVECTORS(128)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::string::Naive<128>>) {
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

TEST_CASE("benchmark vectors symbol() operations - 128 alphabet", "[string_with_rank][!benchmark][128][time][symbol]") {
    auto const& text = generateText<0, 128>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("symbol()")
             .relative(true)
             .batch(text.size());

        call_with_templates<
            ALLRANKVECTORS(128)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::string::Naive<128>>) {
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

TEST_CASE("benchmark vectors rank() operations - 128 alphabet", "[string_with_rank][!benchmark][128][time][rank]") {
    auto const& text = generateText<0, 128>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("rank()")
             .relative(true);

        call_with_templates<
            ALLRANKVECTORS(128)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::string::Naive<128>>) {
                return;
            }

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto rng = ankerl::nanobench::Rng{};

            auto vec = Vector{text};

            bench.run(vector_name, [&]() {
                auto v = vec.rank(rng.bounded(text.size()+1), rng.bounded(128));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        });
    }
}

TEST_CASE("benchmark vectors prefix_rank() operations - 128 alphabet", "[string_with_rank][!benchmark][128][time][prefix_rank]") {
    auto const& text = generateText<0, 128>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("prefix_rank()")
             .relative(true);

        call_with_templates<
            ALLRANKVECTORS(128)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::string::Naive<128>>) {
                return;
            }

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto rng = ankerl::nanobench::Rng{};

            auto vec = Vector{text};

            bench.run(vector_name, [&]() {
                auto v = vec.prefix_rank(rng.bounded(text.size()+1), rng.bounded(128));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        });
    }
}

TEST_CASE("benchmark vectors all_ranks() operations - 128 alphabet", "[string_with_rank][!benchmark][128][time][all_ranks]") {
    auto const& text = generateText<0, 128>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("all_ranks()")
             .relative(true);

        call_with_templates<
            ALLRANKVECTORS(128)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::string::Naive<128>>) {
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

TEST_CASE("benchmark vectors all_ranks_and_prefix_ranks() operations - 128 alphabet", "[string_with_rank][!benchmark][128][time][all_ranks_and_prefix_ranks]") {
    auto const& text = generateText<0, 128>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("all_ranks_and_prefix_ranks()")
             .relative(true);

        call_with_templates<
            ALLRANKVECTORS(128)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::string::Naive<128>>) {
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

TEST_CASE("benchmark vectors in size - alphabet 128", "[string_with_rank][!benchmark][128][size]") {
    auto const& text = generateText<0, 128>();

    SECTION("benchmarking") {
        BenchSize benchSize;
        benchSize.entries[0][2] = "bits/char";
        benchSize.entries[0][3] = "alphabet 128";

        call_with_templates<
            ALLRANKVECTORS(128)>([&]<typename Vector>() {
            if constexpr (std::same_as<Vector, fmindex_collection::string::Naive<128>>) {
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
