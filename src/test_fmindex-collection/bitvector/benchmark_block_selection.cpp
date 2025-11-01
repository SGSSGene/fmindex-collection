// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include <catch2/catch_all.hpp>
#include <cstdlib>
#include <fmindex-collection/utils.h>
#include <fstream>
#include <nanobench.h>

#include "../string/utils.h"
#include "../BenchSize.h"
#include "allBitVectors.h"

//using AllTypes = std::variant<
//    fmc::bitvector::Bitvector2L_512_64k,
//    fmc::bitvector::OptRBBitvector<fmc::bitvector::Bitvector2L_512_64k, fmc::bitvector::Bitvector2L_512_64k>,
//    fmc::bitvector::OptSparseRBBitvector<fmc::bitvector::Bitvector2L_512_64k, fmc::bitvector::Bitvector2L_512_64k>,
//    fmc::bitvector::OptRBBitvector<fmc::bitvector::Bitvector2L_64_64k, fmc::bitvector::Bitvector2L_64_64k>,
//    fmc::bitvector::OptSparseRBBitvector<fmc::bitvector::Bitvector2L_64_64k, fmc::bitvector::Bitvector2L_64_64k>,
//#if defined(FMC_USE_RANKSELECT)
//    RankSelect<5>,
//#endif
//    std::monostate
//>;


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

struct Timer {
    std::chrono::high_resolution_clock::time_point start_;

    Timer() : start_(std::chrono::high_resolution_clock::now()) {}
    void reset() { start_ = std::chrono::high_resolution_clock::now(); }
    double elapsed() const {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(end - start_).count();
    }
};
}


TEST_CASE("benchmark bit vectors for block size selection", "[bitvector][!benchmark][time][rank][block-size]") {
    auto results = std::map<std::string, std::vector<double>>{};
    auto density = 0.10;
    auto blockSize = getEnvOr("BLOCKSIZE", 1);
    auto blockDensity = getEnvOr("BLOCKDENSITY", 0) / 100.;

    auto densitySteps = std::vector<double>{64, 128, 256, 512, 1024, 2048};

        auto text = generateText(density, blockDensity, blockSize);

        using Type1L = std::variant<
            fmc::bitvector::Bitvector1L_64,
            fmc::bitvector::Bitvector1L_128,
            fmc::bitvector::Bitvector1L_256,
            fmc::bitvector::Bitvector1L_512,
            fmc::bitvector::Bitvector1L_1024,
            fmc::bitvector::Bitvector1L_2048,
            std::monostate
        >;
        using Type2L = std::variant<
            fmc::bitvector::Bitvector2L_64_64k,
            fmc::bitvector::Bitvector2L_128_64k,
            fmc::bitvector::Bitvector2L_256_64k,
            fmc::bitvector::Bitvector2L_512_64k,
            fmc::bitvector::Bitvector2L_1024_64k,
            fmc::bitvector::Bitvector2L_2048_64k,
            std::monostate
        >;
        using TypeP1L = std::variant<
            fmc::bitvector::PairedBitvector1L_64,
            fmc::bitvector::PairedBitvector1L_128,
            fmc::bitvector::PairedBitvector1L_256,
            fmc::bitvector::PairedBitvector1L_512,
            fmc::bitvector::PairedBitvector1L_1024,
            fmc::bitvector::PairedBitvector1L_2048,
            std::monostate
        >;
        using TypeP2L = std::variant<
            fmc::bitvector::PairedBitvector2L_64_64k,
            fmc::bitvector::PairedBitvector2L_128_64k,
            fmc::bitvector::PairedBitvector2L_256_64k,
            fmc::bitvector::PairedBitvector2L_512_64k,
            fmc::bitvector::PairedBitvector2L_1024_64k,
            fmc::bitvector::PairedBitvector2L_2048_64k,
            std::monostate
        >;

        using TypeOthers = std::variant<
        #ifdef FMC_USE_PASTA
            FlatRank,
            WideRank,
        #endif
        #ifdef FMC_USE_SDSL
            SDSL_V,
            SDSL_V5,
        #endif
        #ifdef FMC_USE_SUX
            Rank9,
        #endif
        #ifdef FMC_USE_RANKSELECT
            RankSelect<5>,
        #endif
            std::monostate
        >;


        auto measure = [&]<typename AllTypes>(std::string name) {
            auto rng = ankerl::nanobench::Rng{};
            call_with_templates<AllTypes>([&]<typename Vector>() {
                auto vec = Vector{text};

                auto timer = Timer{};
                size_t const numberOfIterations = getEnvOr("ITERATIONS", 1'000);
                for (size_t i{0}; i<numberOfIterations; ++i) {
                    auto v = vec.rank(rng.bounded(text.size()));
                    ankerl::nanobench::doNotOptimizeAway(v);
                };
                auto v = timer.elapsed();
                results[name].push_back(v / (numberOfIterations / 1'000'000'000.));
            });
        };
        measure.operator()<Type1L>("1l");
        measure.operator()<Type2L>("2l");
        measure.operator()<TypeP1L>("p1l");
        measure.operator()<TypeP2L>("p2l");
        measure.operator()<TypeOthers>("others");

    auto nodeTag = std::map<std::string, std::string> {
        {"1l",  "mark=otimes*,mark options={fill=cyan}"},
        {"2l", "mark=otimes*,mark options={fill=olive}"},
        {"p1l", "mark=*,mark options={fill=cyan}"},
        {"p2l", "mark=*,mark options={fill=olive}"},
        {"others", "mark=triangle*,mark options={fill=red}"},
    };
    auto names = std::vector<std::string>{"1l", "2l", "p1l", "p2l", "others"};


    auto ptr = std::getenv("TIKZ_EXPORT");
    if (!ptr) {
        for (auto vector_name : names) {
            auto const& t = results[vector_name];
            for (auto v : t) {
                fmt::print("{: >7.3f} ", v);
            }
            fmt::print("{}\n", vector_name);
        }
    } else {
        for (auto vector_name : names) {
            auto const& t = results[vector_name];
            fmt::print("\\addplot [{}] coordinates {{\n", nodeTag[vector_name]);
            for (size_t i{0}; i < t.size(); ++i) {
                auto x = densitySteps[i];
                auto y = t[i];
                fmt::print("    ({: >3.2f}, {: >7.3f})\n", x, y);
            }
            fmt::print("}};\n");
//            fmt::print("\\addlegendentry{{{}}}\n\n", vector_name);
        }
    }
}
