// SPDX-FileCopyrightText: 2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include <catch2/catch_all.hpp>
#include <cstdlib>
#include <fmindex-collection/utils.h>
#include <fstream>
#include <nanobench.h>

#include "../string/utils.h"
#include "../BenchSize.h"
#include "allBitVectors.h"


namespace {

auto getEnvOr(std::string str, size_t orValue) -> size_t {
    auto ptr = std::getenv(str.c_str());
    if (!ptr) return orValue;
    return std::stoull(ptr);
}
auto generateText(double density, double blockDensity, size_t blockSize, size_t size) -> std::vector<bool> {
    auto rng = ankerl::nanobench::Rng{};

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


TEST_CASE("benchmark bit vectors run time dependent on size", "[bitvector][!benchmark][time][length_vs_time]") {
    auto results = std::map<std::string, std::vector<double>>{};
    auto density = 0.10;
    auto blockSize = getEnvOr("BLOCKSIZE", 1);
    auto blockDensity = getEnvOr("BLOCKDENSITY", 0) / 100.;


    auto sizeSteps = std::vector<size_t>{
        size_t{1}<<24, // 16MB
        size_t{1}<<26, // 64MB
        size_t{1}<<28, // 256MB
        size_t{1}<<30, // 1GB
        size_t{1}<<32, // 4GB
        size_t{1}<<34, // 8GB
    };

    using AllTypes = std::variant<
        fmc::bitvector::Bitvector2L_64_64k,
        fmc::bitvector::Bitvector2L_512_64k,
        fmc::bitvector::PairedBitvector2L_64_64k,
        fmc::bitvector::PairedBitvector2L_512_64k,
        fmc::bitvector::OptSparseRBBitvector<fmc::bitvector::Bitvector2L_64_64k, fmc::bitvector::Bitvector2L_64_64k>,
        fmc::bitvector::OptSparseRBBitvector<fmc::bitvector::Bitvector2L_512_64k, fmc::bitvector::Bitvector2L_512_64k>,
        FlatRank,
        WideRank,
        SDSL_V,
        SDSL_V5,
        Rank9,
        RankSelect<5>,
        std::monostate
    >;

    auto nodeTag = std::map<std::string, std::string> {
        {"2l-64",      "mark=*,mark options={fill=cyan}"},
        {"2l-512",     "mark=otimes*,mark options={fill=cyan}"},
        {"p2l-64",      "mark=*,mark options={fill=cyan}"},
        {"p2l-512",     "mark=otimes*,mark options={fill=cyan}"},
        {"srb-64",     "mark=*,mark options={fill=olive}"},
        {"srb-512",    "mark=otimes*,mark options={fill=olive}"},
        {"FlatRank",   "mark=triangle*,fill=red"},
        {"WideRank",   "mark=triangle*,fill=red!50"},
        {"SDSL_V",     "mark=triangle*,fill=violet"},
        {"SDSL_V5",    "mark=triangle*,mark options={fill=violet!50}"},
        {"Rank9",      "mark=triangle*,mark options={fill=orange}"},
        {"RankSelect", "mark=triangle*,mark options={fill=orange!50}"},
    };
    auto names = std::vector<std::string>{"2l-64", "2l-512", "p2l-64", "p2l-512", "srb-64", "srb-512", "FlatRank", "WideRank", "SDSL_V", "SDSL_V5", "Rank9", "RankSelect"};


    for (auto size : sizeSteps) {
        auto text = generateText(density, blockDensity, blockSize, size);

        auto rng = ankerl::nanobench::Rng{};
        call_with_templates<AllTypes>([&]<typename Vector>() {
            auto vec = Vector{text};
            auto timer = Timer{};
            auto vector_name = getName<Vector>();
            size_t const numberOfIterations = getEnvOr("ITERATIONS", 1'000);
            for (size_t i{0}; i<numberOfIterations; ++i) {
                auto v = vec.rank(rng() % text.size());
                ankerl::nanobench::doNotOptimizeAway(v);
            };
            auto v = timer.elapsed();
            auto index = AllTypes{Vector{}}.index();
            auto name = names[index];
            results[name].push_back(v / (numberOfIterations / 1'000'000'000.));
        });
    }

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
            fmt::print("{: >20}", vector_name);
            auto const& t = results[vector_name];
            for (auto v : t) {
                fmt::print("&  {: >7.3f} ", v);
            }
            fmt::print("\\\\\n");
        }
    }
}
