// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: CC0-1.0
#include <catch2/catch_all.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/optional.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <nanobench.h>
#include <fmindex-collection/bitvector/all.h>
#include <fmindex-collection/suffixarray/SparseArray.h>
#include <fmindex-collection/suffixarray/CompressedSparseArray.h>
#include <map>
#include <unordered_map>

#include "../BenchSize.h"


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

//    using Load = std::array<size_t, 8>;
    using Load = uint64_t;

    auto rng = ankerl::nanobench::Rng{};
    static auto indicatorBV = fmc::bitvector::PairedL0L1_512_64kBitvector{std::views::iota(size_t{}, size) | std::views::transform([&](size_t) -> bool {
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
                    buffer[i] = Load{i};
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
                    buffer[i] = Load{i};
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
            auto sparseArray = fmc::suffixarray::SparseArray<Load>{
                std::views::iota(size_t{0}, indicatorBV.size()) | std::views::transform([&](size_t i) -> std::optional<Load> {
                    if (!indicatorBV.symbol(i)) return std::nullopt;
                    return i;
                })
            };
            bench.run("sparse array", [&]() {
                auto pos = rng.bounded(size);
                auto v = sparseArray.value(pos);
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        }

        if ((mode & 0x10) != 0) {
            auto sparseArray = fmc::suffixarray::CompressedSparseArray {
                std::views::iota(size_t{0}, indicatorBV.size()) | std::views::transform([&](size_t i) -> std::optional<Load> {
                    if (!indicatorBV.symbol(i)) return std::nullopt;
                    return i;
                })
            };
            bench.run("compressed sparse array", [&]() {
                auto pos = rng.bounded(size);
                auto v = sparseArray.value(pos);
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        }

        if ((mode & 0x20) != 0) {
            auto sparseArray = fmc::suffixarray::CompressedSparseArrayV2 {
                std::views::iota(size_t{0}, indicatorBV.size()) | std::views::transform([&](size_t i) -> std::optional<Load> {
                    if (!indicatorBV.symbol(i)) return std::nullopt;
                    return i;
                })
            };
            bench.run("compressed sparse array v2", [&]() {
                auto pos = rng.bounded(size);
                auto v = sparseArray.value(pos);
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        }

        if ((mode & 0x40) != 0) {
            auto sparseArray = fmc::suffixarray::CompressedSparseArrayV3 {
                std::views::iota(size_t{0}, indicatorBV.size()) | std::views::transform([&](size_t i) -> std::optional<Load> {
                    if (!indicatorBV.symbol(i)) return std::nullopt;
                    return i;
                })
            };
            bench.run("compressed sparse array v3", [&]() {
                auto pos = rng.bounded(size);
                auto v = sparseArray.value(pos);
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        }
        if ((mode & 0x80) != 0) {
            auto sparseArray = fmc::suffixarray::CompressedSparseArrayV4 {
                std::views::iota(size_t{0}, indicatorBV.size()) | std::views::transform([&](size_t i) -> std::optional<Load> {
                    if (!indicatorBV.symbol(i)) return std::nullopt;
                    return i;
                })
            };
            bench.run("compressed sparse array v4", [&]() {
                auto pos = rng.bounded(size);
                auto v = sparseArray.value(pos);
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        }

        if ((mode & 0x100) != 0) {
            auto sparseArray = fmc::suffixarray::CompressedSparseArrayV5 {
                std::views::iota(size_t{0}, indicatorBV.size()) | std::views::transform([&](size_t i) -> std::optional<Load> {
                    if (!indicatorBV.symbol(i)) return std::nullopt;
                    return i;
                })
            };
            bench.run("compressed sparse array v5", [&]() {
                auto pos = rng.bounded(size);
                auto v = sparseArray.value(pos);
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        }

    }
}

TEST_CASE("benchmark sparse array - size", "[sparsearray][!benchmark][size]") {
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

    BenchSize benchSize;
    benchSize.baseSize = 64 / rate;
    benchSize.entries[0][2] = "bits/entry";
    benchSize.entries[0][4] = "compressed sparse array";



//    using Load = std::array<size_t, 8>;
    using Load = uint64_t;

    auto rng = ankerl::nanobench::Rng{};
    static auto indicatorBV = fmc::bitvector::PairedL0L1_512_64kBitvector{std::views::iota(size_t{}, size) | std::views::transform([&](size_t) -> bool {
        return rng.bounded(rate) == 0;
    })};

//    auto nbrOfOnes = indicatorBV.rank(indicatorBV.size());

    auto mode = []() -> size_t {
        auto ptr = std::getenv("SPARSEARRAYMODE");
        if (ptr) {
            return std::stoull(ptr);
        }
        return 0;
    }();

    SECTION("benchmarking") {
        if ((mode & 0x01) != 0) {
            auto buffer = std::unordered_map<uint64_t, Load>{};
            for (size_t i{}; i < indicatorBV.size(); ++i) {
                if (indicatorBV.symbol(i)) {
                    buffer[i] = Load{i};
                }
            }

/*            auto size = [&]() {
                auto ofs     = std::stringstream{};
                auto archive = cereal::BinaryOutputArchive{ofs};
                archive(buffer);
                return ofs.str().size();
            }();

            benchSize.addEntry({
                .name = "unordered_map",
                .size = size,
                .text_size = nbrOfOnes*64/rate,
                .bits_per_char = (size*8)/double(nbrOfOnes)
            });*/
        }

        if ((mode & 0x02) != 0) {
            auto buffer = std::map<uint64_t, Load>{};
            for (size_t i{}; i < indicatorBV.size(); ++i) {
                if (indicatorBV.symbol(i)) {
                    buffer[i] = Load{i};
                }
            }
/*            auto size = [&]() {
                auto ofs     = std::stringstream{};
                auto archive = cereal::BinaryOutputArchive{ofs};
                archive(buffer);
                return ofs.str().size();
            }();

            benchSize.addEntry({
                .name = "map",
                .size = size,
                .text_size = nbrOfOnes,
                .bits_per_char = (size*8)/double(nbrOfOnes)
            });*/
        }

        if ((mode & 0x04) != 0) {
            auto buffer = std::vector<std::optional<Load>>{};
            buffer.resize(indicatorBV.size());
            for (size_t i{}; i < indicatorBV.size(); ++i) {
                if (indicatorBV.symbol(i)) {
                    buffer[i] = Load{};
                }
            }

/*            auto size = [&]() {
                auto ofs     = std::stringstream{};
                auto archive = cereal::BinaryOutputArchive{ofs};
                archive(buffer);
                return ofs.str().size();
            }();

            benchSize.addEntry({
                .name = "vector",
                .size = size,
                .text_size = nbrOfOnes*64/rate,
                .bits_per_char = (size*8)/double(nbrOfOnes)
            });*/
        }

        if ((mode & 0x08) != 0) {
            auto sparseArray = fmc::suffixarray::SparseArray<Load>{
                std::views::iota(size_t{0}, indicatorBV.size()) | std::views::transform([&](size_t i) -> std::optional<Load> {
                    if (!indicatorBV.symbol(i)) return std::nullopt;
                    return i;
                })
            };
            (void)sparseArray;
/*            auto size = [&]() {
                auto ofs     = std::stringstream{};
                auto archive = cereal::BinaryOutputArchive{ofs};
                archive(sparseArray);
                return ofs.str().size();
            }();

            benchSize.addEntry({
                .name = "sparse array",
                .size = size,
                .text_size = nbrOfOnes,
                .bits_per_char = (size*8)/double(nbrOfOnes)
            });*/
        }

        if ((mode & 0x10) != 0) {
            auto sparseArray = fmc::suffixarray::CompressedSparseArray {
                std::views::iota(size_t{0}, indicatorBV.size()) | std::views::transform([&](size_t i) -> std::optional<Load> {
                    if (!indicatorBV.symbol(i)) return std::nullopt;
                    return i;
                })
            };
            (void)sparseArray;
/*            auto size = [&]() {
                auto ofs     = std::stringstream{};
                auto archive = cereal::BinaryOutputArchive{ofs};
                archive(sparseArray);
                return ofs.str().size();
            }();

            benchSize.addEntry({
                .name = "compressed sparse array",
                .size = size,
                .text_size = nbrOfOnes,
                .bits_per_char = (size*8)/double(nbrOfOnes)
            });*/
        }

        if ((mode & 0x20) != 0) {
            auto sparseArray = fmc::suffixarray::CompressedSparseArrayV2 {
                std::views::iota(size_t{0}, indicatorBV.size()) | std::views::transform([&](size_t i) -> std::optional<Load> {
                    if (!indicatorBV.symbol(i)) return std::nullopt;
                    return i;
                })
            };
            (void)sparseArray;
/*            auto size = [&]() {
                auto ofs     = std::stringstream{};
                auto archive = cereal::BinaryOutputArchive{ofs};
                archive(sparseArray);
                return ofs.str().size();
            }();

            benchSize.addEntry({
                .name = "compressed sparse array v2",
                .size = size,
                .text_size = nbrOfOnes,
                .bits_per_char = (size*8)/double(nbrOfOnes)
            });*/
        }
        if ((mode & 0x40) != 0) {
            auto sparseArray = fmc::suffixarray::CompressedSparseArrayV3 {
                std::views::iota(size_t{0}, indicatorBV.size()) | std::views::transform([&](size_t i) -> std::optional<Load> {
                    if (!indicatorBV.symbol(i)) return std::nullopt;
                    return i;
                })
            };
            (void)sparseArray;
/*            auto size = [&]() {
                auto ofs     = std::stringstream{};
                auto archive = cereal::BinaryOutputArchive{ofs};
                archive(sparseArray);
                return ofs.str().size();
            }();

            benchSize.addEntry({
                .name = "compressed sparse array v3",
                .size = size,
                .text_size = nbrOfOnes,
                .bits_per_char = (size*8)/double(nbrOfOnes)
            });*/
        }
        if ((mode & 0x80) != 0) {
            auto sparseArray = fmc::suffixarray::CompressedSparseArrayV4 {
                std::views::iota(size_t{0}, indicatorBV.size()) | std::views::transform([&](size_t i) -> std::optional<Load> {
                    if (!indicatorBV.symbol(i)) return std::nullopt;
                    return i;
                })
            };
            (void)sparseArray;
/*            auto size = [&]() {
                auto ofs     = std::stringstream{};
                auto archive = cereal::BinaryOutputArchive{ofs};
                archive(sparseArray);
                return ofs.str().size();
            }();

            benchSize.addEntry({
                .name = "compressed sparse array v3",
                .size = size,
                .text_size = nbrOfOnes,
                .bits_per_char = (size*8)/double(nbrOfOnes)
            });*/
        }

        if ((mode & 0x100) != 0) {
            auto sparseArray = fmc::suffixarray::CompressedSparseArrayV5 {
                std::views::iota(size_t{0}, indicatorBV.size()) | std::views::transform([&](size_t i) -> std::optional<Load> {
                    if (!indicatorBV.symbol(i)) return std::nullopt;
                    return i;
                })
            };
            (void)sparseArray;
/*            auto size = [&]() {
                auto ofs     = std::stringstream{};
                auto archive = cereal::BinaryOutputArchive{ofs};
                archive(sparseArray);
                return ofs.str().size();
            }();

            benchSize.addEntry({
                .name = "compressed sparse array v3",
                .size = size,
                .text_size = nbrOfOnes,
                .bits_per_char = (size*8)/double(nbrOfOnes)
            });*/
        }

    }
}
