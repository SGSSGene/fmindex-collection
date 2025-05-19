// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: CC0-1.0

#include "../string/utils.h"

TEST_CASE("benchmark binary search on large numbers", "[!benchmark][binary-search][time]") {
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
            auto e = rng.bounded(data.back());
            auto v = std::ranges::lower_bound(data, e);
            ankerl::nanobench::doNotOptimizeAway(*v);
        });
    }
}

