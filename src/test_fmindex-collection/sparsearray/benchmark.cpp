// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: CC0-1.0
#include <catch2/catch_all.hpp>
#include <nanobench.h>
#include <fmindex-collection/bitvector/all.h>
#include <map>

TEST_CASE("benchmark sparse array", "[sparsearray][!benchmark][time]") {
    auto size = []() -> size_t {
        auto ptr = std::getenv("SPARSEARRAYSIZE");
        if (ptr) {
            return std::stoull(ptr);
        }
        #ifdef NDEBUG
            return 10'000'000;
        #else
            return 100'000;
        #endif
    }();

    auto rate = []() -> size_t {
        auto ptr = std::getenv("SPARSEARRAYRATE");
        if (ptr) {
            return std::stoull(ptr);
        }
        return 16;
    }();

    using Load = std::array<size_t, 8>;

    auto rng = ankerl::nanobench::Rng{};
    static auto indicatorBV = fmindex_collection::bitvector::PairedL0L1_512_64kBitvector{std::views::iota(size_t{}, size) | std::views::transform([&](size_t) -> bool {
        return rng.bounded(rate) == 0;
    })};

    auto mode = []() -> size_t {
        auto ptr = std::getenv("SPARSEARRAYMODE");
        if (ptr) {
            return std::stoull(ptr);
        }
        return 0;
    }();

    auto bench = ankerl::nanobench::Bench{};
    bench.title("access()")
         .relative(true);
    SECTION("benchmarking") {

        if ((mode & 0x01) != 0) {
            auto buffer = std::unordered_map<uint64_t, Load>{};
            for (size_t i{}; i < indicatorBV.size(); ++i) {
                if (indicatorBV.symbol(i)) {
                    buffer[i] = Load{};
                }
            }
            bench.run("unordered_map", [&]() {
                auto pos = rng.bounded(size);
                auto iter = buffer.find(pos);
                if (iter != buffer.end()) {
                    auto v = iter->second;
                    ankerl::nanobench::doNotOptimizeAway(v);
                }
            });
        }

        if ((mode & 0x02) != 0) {
            auto buffer = std::map<uint64_t, Load>{};
            for (size_t i{}; i < indicatorBV.size(); ++i) {
                if (indicatorBV.symbol(i)) {
                    buffer[i] = Load{};
                }
            }
            bench.run("map", [&]() {
                auto pos = rng.bounded(size);
                auto iter = buffer.find(pos);
                if (iter != buffer.end()) {
                    auto v = iter->second;
                    ankerl::nanobench::doNotOptimizeAway(v);
                }
            });
        }

        if ((mode & 0x04) != 0) {
            auto buffer = std::vector<std::optional<Load>>{};
            buffer.resize(indicatorBV.size());
            for (size_t i{}; i < indicatorBV.size(); ++i) {
                if (indicatorBV.symbol(i)) {
                    buffer[i] = Load{};
                }
            }
            bench.run("vector", [&]() {
                auto pos = rng.bounded(size);
                if (buffer[pos]) {
                    auto v = *buffer[pos];
                    ankerl::nanobench::doNotOptimizeAway(v);
                }
            });
        }

        if ((mode & 0x08) != 0) {
            auto bv = indicatorBV;
            auto values = std::vector<Load>{};
            values.reserve(indicatorBV.rank(indicatorBV.size()));
            for (size_t i{0}; i < size; ++i) {
                if (bv.symbol(i)) {
                    values.emplace_back();
                }
            }
            bench.run("sparse array", [&]() {
                auto pos = rng.bounded(size);
                auto b = bv.symbol(pos);
                if (b) {
                    auto r = bv.rank(pos);
                    auto v = values[r];
                    ankerl::nanobench::doNotOptimizeAway(v);
                }
            });
        }
    }
}
