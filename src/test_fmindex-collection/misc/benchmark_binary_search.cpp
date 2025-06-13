// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: CC0-1.0

#include "../string/utils.h"

TEST_CASE("benchmark binary search on large numbers", "[!benchmark][binary-search][uint64_t][time]") {
    auto size = []() -> size_t {
        auto ptr = std::getenv("STRINGSIZE");
        #ifdef NDEBUG
        size_t defVal = 1'000'000;
        #else
        size_t defVal = 1'000;
        #endif

        return ankerl::nanobench::detail::strToUInt(ptr, defVal);
    }();


    static auto data = [&]() -> std::vector<uint64_t> {
        auto data = std::vector<uint64_t>{};
        data.resize(size);
        auto rng = ankerl::nanobench::Rng{};
        for (size_t i{0}; i < size; ++i) {
            data[i] = rng();
        }
        std::ranges::sort(data);
        return data;
    }();



    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("lower_bound")
             .relative(true);
        auto rng = ankerl::nanobench::Rng{};
        bench.run("lower_bound 64bit", [&]() {
            auto e = rng();
            auto v = std::ranges::lower_bound(data, e);
            ankerl::nanobench::doNotOptimizeAway(*v);
        });
    }
}

TEST_CASE("benchmark binary search on large numbers", "[!benchmark][binary-search][uint32_t][time]") {
    auto size = []() -> size_t {
        auto ptr = std::getenv("STRINGSIZE");
        #ifdef NDEBUG
        size_t defVal = 1'000'000;
        #else
        size_t defVal = 1'000;
        #endif

        return ankerl::nanobench::detail::strToUInt(ptr, defVal);
    }();


    static auto data = [&]() -> std::vector<uint32_t> {
        auto data = std::vector<uint32_t>{};
        data.resize(size);
        auto rng = ankerl::nanobench::Rng{};
        for (size_t i{0}; i < size; ++i) {
            data[i] = rng.bounded(std::numeric_limits<uint32_t>::max());
        }
        std::ranges::sort(data);
        return data;
    }();



    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("lower_bound")
             .relative(true);
        auto rng = ankerl::nanobench::Rng{};
        bench.run("lower_bound 32bit", [&]() {
            auto e = rng.bounded(std::numeric_limits<uint32_t>::max());
            auto v = std::ranges::lower_bound(data, e);
            ankerl::nanobench::doNotOptimizeAway(*v);
        });
    }
}

TEST_CASE("benchmark binary search on large numbers", "[!benchmark][binary-search][double][time]") {
    auto size = []() -> size_t {
        auto ptr = std::getenv("STRINGSIZE");
        #ifdef NDEBUG
        size_t defVal = 1'000'000;
        #else
        size_t defVal = 1'000;
        #endif

        return ankerl::nanobench::detail::strToUInt(ptr, defVal);
    }();


    static auto data = [&]() -> std::vector<double> {
        auto data = std::vector<double>{};
        data.resize(size);
        auto rng = ankerl::nanobench::Rng{};
        for (size_t i{0}; i < size; ++i) {
            data[i] = rng.uniform01();
        }
        std::ranges::sort(data);
        return data;
    }();



    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("lower_bound")
             .relative(true);
        auto rng = ankerl::nanobench::Rng{};
        bench.run("lower_bound double", [&]() {
            auto e = rng.uniform01();
            auto v = std::ranges::lower_bound(data, e);
            ankerl::nanobench::doNotOptimizeAway(*v);
        });
    }
}

TEST_CASE("benchmark binary search on large numbers", "[!benchmark][binary-search][float][time]") {
    auto size = []() -> size_t {
        auto ptr = std::getenv("STRINGSIZE");
        #ifdef NDEBUG
        size_t defVal = 1'000'000;
        #else
        size_t defVal = 1'000;
        #endif

        return ankerl::nanobench::detail::strToUInt(ptr, defVal);
    }();


    static auto data = [&]() -> std::vector<float> {
        auto data = std::vector<float>{};
        data.resize(size);
        auto rng = ankerl::nanobench::Rng{};
        for (size_t i{0}; i < size; ++i) {
            data[i] = rng.uniform01();
        }
        std::ranges::sort(data);
        return data;
    }();



    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("lower_bound")
             .relative(true);
        auto rng = ankerl::nanobench::Rng{};
        bench.run("lower_bound float", [&]() {
            auto e = rng.uniform01();
            auto v = std::ranges::lower_bound(data, e);
            ankerl::nanobench::doNotOptimizeAway(*v);
        });
    }
}



TEST_CASE("benchmark binary search on large numbers and tuples", "[!benchmark][binary-search][tuple][time]") {
    auto size = []() -> size_t {
        auto ptr = std::getenv("STRINGSIZE");
        #ifdef NDEBUG
        size_t defVal = 1'000'000;
        #else
        size_t defVal = 1'000;
        #endif

        return ankerl::nanobench::detail::strToUInt(ptr, defVal);
    }();


    static auto data = [&]() -> std::vector<std::tuple<uint64_t, uint64_t>> {
        auto data = std::vector<std::tuple<uint64_t, uint64_t>>{};
        data.resize(size);
        auto rng = ankerl::nanobench::Rng{};
        for (size_t i{0}; i < size; ++i) {
            data[i] = std::make_tuple(rng(), rng());
        }
        std::ranges::sort(data);
        return data;
    }();



    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("lower_bound")
             .relative(true);
        auto rng = ankerl::nanobench::Rng{};
        bench.run("lower_bound tuple", [&]() {
            auto e = std::make_tuple(rng(), rng());
            auto v = std::ranges::lower_bound(data, e);
            ankerl::nanobench::doNotOptimizeAway(*v);
        });
    }
}

