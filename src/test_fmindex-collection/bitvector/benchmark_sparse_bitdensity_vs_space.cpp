// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: CC0-1.0
#include <catch2/catch_all.hpp>
#include <cstdlib>
#include <fmindex-collection/utils.h>
#include <fstream>
#include <nanobench.h>

#include "../string/utils.h"
#include "../BenchSize.h"
#include "allBitVectors.h"

using AllTypes = std::variant<
    fmc::bitvector::Bitvector2L_512_64k,
    fmc::bitvector::OptRBBitvector<fmc::bitvector::Bitvector2L_512_64k, fmc::bitvector::Bitvector2L_512_64k>,
    fmc::bitvector::OptSparseRBBitvector<fmc::bitvector::Bitvector2L_512_64k, fmc::bitvector::Bitvector2L_512_64k>,
    fmc::bitvector::OptRBBitvector<fmc::bitvector::Bitvector2L_64_64k, fmc::bitvector::Bitvector2L_64_64k>,
    fmc::bitvector::OptSparseRBBitvector<fmc::bitvector::Bitvector2L_64_64k, fmc::bitvector::Bitvector2L_64_64k>,
#if defined(FMC_USE_RANKSELECT)
    RankSelect<5>,
#endif
    std::monostate
>;

namespace {
constexpr auto nodeTag = std::array<std::string_view, std::variant_size_v<AllTypes>> {
    "mark=triangle*,mark options={fill=red}",
    "mark=otimes*,mark options={fill=cyan}",
    "mark=otimes*,mark options={fill=olive}",
    "mark=*,mark options={fill=cyan}",
    "mark=*,mark options={fill=olive}",
#if defined(FMC_USE_RANKSELECT)
    "mark=triangle*,mark options={fill=violet}",
#endif
    ""
};
}

namespace {

auto getEnvOr(std::string str, size_t orValue) -> size_t {
    auto ptr = std::getenv(str.c_str());
    if (!ptr) return orValue;
    return std::stoull(ptr);
}
/*auto generateText(double density) -> std::vector<bool> {
    auto rng = ankerl::nanobench::Rng{};

    auto size = getEnvOr("BITVECTORSIZE", 10'000'000);

    auto text = std::vector<bool>{};
    for (size_t i{0}; i<size; ++i) {
        text.push_back(rng.bounded(100'000) < density*100'000.);
    }
    return text;
}*/
auto generateText(double density, double blockDensity, size_t blockSize) -> std::vector<bool> {
    auto rng = ankerl::nanobench::Rng{};

    auto size = getEnvOr("BITVECTORSIZE", 10'000'000);

    auto text = std::vector<bool>{};
    for (size_t i{0}; i+blockSize-1 < size; i += blockSize) {
        auto isBlock = (rng.bounded(100'000) < blockDensity*100'000);
        if (isBlock) {
            auto v = rng.bounded(100'000) < density*100'000;
            for (size_t j{0}; j < blockSize; ++j) {
                text.push_back(v);
            }
        } else {
            for (size_t j{0}; j < blockSize; ++j) {
                text.push_back(rng.bounded(100'000) < density*100'000);
            }
        }

    }
    while (text.size() < size) {
        text.push_back(rng.bounded(100'000) < density*100'000);
    }
    return text;
}

}

TEST_CASE("benchmark bitdensity against space", "[sparse-bitvector][!benchmark][bitdensity_vs_space]") {
    auto results = std::map<std::string, std::vector<double>>{};
    auto densitySteps = std::vector<double>{0.01, 0.02, 0.03, 0.04, 0.05, 0.10, 0.15, 0.20, 0.25, 0.3, 0.35, 0.40, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99};

    auto blockSize = getEnvOr("BLOCKSIZE", 1);
    auto blockDensity = getEnvOr("BLOCKDENSITY", 0) / 100.;
    for (auto density : densitySteps) {
        auto text = generateText(density, blockDensity, blockSize);
        call_with_templates<AllTypes>([&]<typename Vector>() {
            auto vector_name = getName<Vector>();
            auto rng = ankerl::nanobench::Rng{};
            auto vec = Vector{text};
            auto spaceUsage = [&]() {
                if constexpr (requires() { vec.space_usage(); }) {
                    auto s = vec.space_usage();
                    return s;
                } else {
                    auto ofs     = std::stringstream{};
                    auto archive = cereal::BinaryOutputArchive{ofs};
                    archive(vec);
                    auto s = ofs.str().size();
                    return s;
                }
            }();
            auto bits_per_bit = (spaceUsage*8) / double(text.size());
            results[vector_name].push_back(bits_per_bit);
        });
    }

    auto ptr = std::getenv("TIKZ_EXPORT");
    if (!ptr) {
        call_with_templates<AllTypes>([&]<typename Vector>() {
            auto vector_name = getName<Vector>();
            auto const& t = results[vector_name];
            for (auto v : t) {
                fmt::print("{: >7.3f} ", v);
            }
            fmt::print("{}\n", vector_name);
        });
    } else {
        call_with_templates<AllTypes>([&]<typename Vector>() {
            auto vector_name = getName<Vector>();
            auto const& t = results[vector_name];
            auto index = AllTypes{Vector{}}.index();
            fmt::print("\\addplot [{}] coordinates {{\n", nodeTag[index]);
            for (size_t i{0}; i < t.size(); ++i) {
                auto x = densitySteps[i];
                auto y = t[i];
                fmt::print("    ({: >3.2f}, {: >7.3f})\n", x, y);
            }
            fmt::print("}};\n");
//            fmt::print("\\addlegendentry{{{}}}\n\n", vector_name);
        });
    }
}
