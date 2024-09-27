// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include <catch2/catch_all.hpp>
#include <fmindex-collection/utils.h>
#include <fstream>
#include <nanobench.h>

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
    Bench bench_symbol{"symbol()"};
    Bench bench_ctor{"c'tor"};
};
}

static auto benchs = Benchs{};
TEMPLATE_TEST_CASE("benchmark bit vectors ctor run times", "[BitVector][!benchmark][time][ctor][.]", ALLBITVECTORS) {
    using Vector = TestType;
    auto vector_name = getName<Vector>();
    INFO(vector_name);

    auto& [bench_rank, bench_symbol, bench_ctor] = benchs;

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

        size_t minEpochIterations = 2'000'000;
        minEpochIterations = 1;

        bench_ctor.minEpochIterations(minEpochIterations).batch(text.size()).run(vector_name, [&]() {
            auto vec = Vector{text};
            ankerl::nanobench::doNotOptimizeAway(vec);
        });
    }
}

static auto seed = ankerl::nanobench::Rng{0}();
static auto text_time = []() -> std::vector<uint8_t> {
    auto rng = ankerl::nanobench::Rng{seed};

    auto text = std::vector<uint8_t>{};
    #ifdef NDEBUG
    for (size_t i{0}; i<10'000'000; ++i) {
    #else
    for (size_t i{0}; i<100'000; ++i) {
    #endif
        text.push_back(rng.bounded(4) == 0);
    }
    return text;
}();

TEMPLATE_TEST_CASE("benchmark bit vectors rank and symbol run times", "[BitVector][!benchmark][time]", ALLBITVECTORS) {
    using Vector = TestType;
    auto vector_name = getName<Vector>();
    INFO(vector_name);

    auto& [bench_rank, bench_symbol, bench_ctor] = benchs;

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        auto rng = ankerl::nanobench::Rng{};

        auto const& text = text_time;
        auto vec = Vector{text};

        size_t minEpochIterations = 10'000'000;
//        minEpochIterations = 1;

        bench_symbol.minEpochIterations(minEpochIterations).run(vector_name, [&]() {
            auto v = vec.symbol(rng.bounded(text.size()));
            ankerl::nanobench::doNotOptimizeAway(v);
        });
        bench_rank.minEpochIterations(minEpochIterations).run(vector_name, [&]() {
            auto v = vec.rank(rng.bounded(text.size()));
            ankerl::nanobench::doNotOptimizeAway(v);
        });
    }
}

namespace {
struct BenchSize {
    struct Entry {
        std::string name;
        size_t      size;
        size_t      text_size;
        double      bits_per_char;
        double      relative{};
    };

    void addEntry(Entry e) {
        entries.push_back(e);
        entries.back().relative = double(entries.back().size) / entries[0].size * 100.;
    }

    std::vector<Entry> entries;

    ~BenchSize() {
        if (entries.empty()) return;
        auto lines = std::vector<std::string>{};
        lines.resize(entries.size());

        auto addColumn = [&]<typename T>(std::format_string<T&, size_t&> fmt, T BenchSize::Entry::*ptr) {
            size_t longestLine{};
            for (size_t i{0}; i < lines.size(); ++i) {
                size_t x = 1;
                auto t = std::format(fmt, entries[i].*ptr, x);
                longestLine = std::max(longestLine, t.size());
            }
            for (size_t i{0}; i < lines.size(); ++i) {
                auto t = std::format(fmt, entries[i].*ptr, longestLine);
                lines[i] = std::format("{} | {}", lines[i], t);
            }
        };

        addColumn("{:> {}.1f}%", &BenchSize::Entry::relative);
        addColumn("{:> {}}", &BenchSize::Entry::size);
        addColumn("{:> {}.3f}", &BenchSize::Entry::bits_per_char);
        addColumn("{:<{}}", &BenchSize::Entry::name);

        for (size_t i{0}; i < lines.size(); ++i) {
            std::cout << lines[i] << "\n";
        }
    }
};
BenchSize benchSize;
}

TEMPLATE_TEST_CASE("benchmark bit vectors memory consumption", "[BitVector][!benchmark][size]", ALLBITVECTORS) {
    using Vector = TestType;
    auto vector_name = getName<Vector>();
    INFO(vector_name);

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

    }
}
