// SPDX-FileCopyrightText: 2026 Simon Gene Gottlieb
// SPDX-License-Identifier: CC0-1.0

#include "../string/allStrings.h"
#include "../string/utils.h"
#include <fmindex-collection/fmindex/LinearFMIndexCursor.h>

TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][uint64]") {

    auto data = std::vector<uint64_t>{};
    auto rng = ankerl::nanobench::Rng{};
    auto rng64 = [&]() {
        return uint64_t{rng()} << 32 | rng();
    };
    {
        // generates string with values between 1-4
        auto size = []() -> size_t {
            auto ptr = std::getenv("STRINGSIZE");
            if (ptr) {
                return std::stoull(ptr);
            }
            #ifdef NDEBUG
                return 1'000'000;
            #else
                return 1'000;
            #endif
        }();

        for (size_t i{0}; i<size; ++i) {
            data.push_back(rng64());
        }
    }
    std::sort(data.begin(), data.end());

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench
            .title("search uint64_t")
            .relative(true);

        bench.run("std::equal_range", [&]() {
            auto [iter1, iter2] = std::equal_range(data.begin(), data.end(), rng64());

            ankerl::nanobench::doNotOptimizeAway(iter1);
            ankerl::nanobench::doNotOptimizeAway(iter2);
        });

        auto const& sequences = reinterpret_cast<std::vector<std::array<uint8_t, 8>> const&>(data);

        auto index = fmc::LinearFMIndex<256, fmc::string::PairedFlattenedBitvectors_512_64k>{sequences};
        bench.run("LinearFMIndex", [&]() {
            auto cur = fmc::LinearFMIndexCursor{index};
            auto v = rng64();
            auto const& vs = reinterpret_cast<std::array<uint8_t, 8> const&>(v);
            for (size_t i{0}; i < 8; ++i) {
                cur = cur.extendLeft(vs[7-i]);
            }
            ankerl::nanobench::doNotOptimizeAway(cur);
        });

    }
}


template <size_t Bytes, size_t Sigma>
static void benchmark() {
    using Value = std::array<uint8_t, Bytes>;
    auto data = std::vector<Value>{};
    auto rng = ankerl::nanobench::Rng{};
    auto rngV = [&]() {
        auto value = Value{};
        for (auto& v : value) {
            v = rng.bounded(Sigma);
        }
        return value;
    };
    {
        // generates string with values between 1-4
        auto size = []() -> size_t {
            auto ptr = std::getenv("STRINGSIZE");
            if (ptr) {
                return std::stoull(ptr);
            }
            #ifdef NDEBUG
                return 1'000'000;
            #else
                return 1'000;
            #endif
        }();

        for (size_t i{0}; i<size; ++i) {
            data.push_back(rngV());
        }
    }
    auto bench = ankerl::nanobench::Bench{};
    bench
        .title(std::string{"search "} + std::to_string(Bytes) + "b  sigma="+std::to_string(Sigma))
        .relative(true);

    {
        auto index = fmc::LinearFMIndex<Sigma, fmc::string::PairedFlattenedBitvectors_512_64k>{data};
        bench.run("LinearFMIndex", [&]() {
            auto cur = fmc::LinearFMIndexCursor{index};
            auto vs = rngV();
            for (size_t i{0}; i < std::tuple_size_v<Value>; ++i) {
                cur = cur.extendLeft(vs[std::tuple_size_v<Value>-1-i]);
            }
            ankerl::nanobench::doNotOptimizeAway(cur);
        });
    }
    {
        auto indices = std::vector<uint64_t>{};
        indices.resize(data.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), [&](auto lhs, auto rhs) {
            return data[lhs] < data[rhs];
        });
        bench.run("std::equal_range - indirect", [&]() {
            auto v = rngV();
            auto iter1 = std::lower_bound(indices.begin(), indices.end(), v, [&](auto const& lhs, auto const& rhs) {
                return data[lhs] < rhs;
            });
            auto iter2 = std::upper_bound(indices.begin(), indices.end(), v, [&](auto const& lhs, auto const& rhs) {
                return lhs < data[rhs];
            });


            ankerl::nanobench::doNotOptimizeAway(iter1);
            ankerl::nanobench::doNotOptimizeAway(iter2);
        });
    }
    {
        std::sort(data.begin(), data.end());
        bench.run("std::equal_range", [&]() {
            auto [iter1, iter2] = std::equal_range(data.begin(), data.end(), rngV());

            ankerl::nanobench::doNotOptimizeAway(iter1);
            ankerl::nanobench::doNotOptimizeAway(iter2);
        });
    }
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][1b][s4]") {
    benchmark<1, 4>();
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][2b][s4]") {
    benchmark<2, 4>();
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][4b][s4]") {
    benchmark<4, 4>();
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][8b][s4]") {
    benchmark<8, 4>();
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][16b][s4]") {
    benchmark<16, 4>();
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][32b][s4]") {
    benchmark<32, 4>();
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][64b][s4]") {
    benchmark<64, 4>();
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][128b][s4]") {
    benchmark<128, 4>();
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][1b][s16]") {
    benchmark<1, 16>();
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][2b][s16]") {
    benchmark<2, 16>();
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][4b][s16]") {
    benchmark<4, 16>();
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][8b][s16]") {
    benchmark<8, 16>();
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][16b][s16]") {
    benchmark<16, 16>();
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][32b][s16]") {
    benchmark<32, 16>();
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][64b][s16]") {
    benchmark<64, 16>();
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][128b][s16]") {
    benchmark<128, 16>();
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][1b][s256]") {
    benchmark<1, 256>();
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][2b][s256]") {
    benchmark<2, 256>();
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][4b][s256]") {
    benchmark<4, 256>();
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][8b][s256]") {
    benchmark<8, 256>();
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][16b][s256]") {
    benchmark<16, 256>();
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][32b][s256]") {
    benchmark<32, 256>();
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][64b][s256]") {
    benchmark<64, 256>();
}
TEST_CASE("benchmark search std::equal_range", "[linearfmindex][!benchmark][128b][s256]") {
    benchmark<128, 256>();
}
