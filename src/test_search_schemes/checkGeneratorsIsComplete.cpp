// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include <catch2/catch_all.hpp>
#include <search_schemes/generator/backtracking.h>
#include <search_schemes/generator/bestKnown.h>
#include <search_schemes/generator/h2.h>
#include <search_schemes/generator/kianfar.h>
#include <search_schemes/generator/kucherov.h>
#include <search_schemes/generator/optimum.h>
#include <search_schemes/generator/pigeon.h>
#include <search_schemes/generator/suffixFilter.h>
#include <search_schemes/generator/zeroOnesZero.h>
#include <search_schemes/isComplete.h>
#include <search_schemes/isValid.h>

namespace ss = search_schemes;
namespace gen = ss::generator;

TEST_CASE("check search scheme generator backtracking for completness", "[isComplete][backtracking]") {
    for (size_t N{1}; N < 10; ++N) { // Number of pieces
        INFO("N " << N);
        for (size_t minK{0}; minK < 5ull; ++minK) {
            INFO("minK " << minK);
            for (size_t maxK{minK}; maxK < 5ull; ++maxK) {
                INFO("maxK " << maxK);
                CHECK(ss::isComplete(gen::backtracking(N, minK, maxK), minK, maxK));
            }
        }
    }
}

TEST_CASE("check search scheme generator bestKnown for completness", "[isComplete][bestKnown]") {
    for (size_t N{1}; N < 10; ++N) { // Number of pieces
        INFO("N " << N);
        for (size_t minK{0}; minK < 5ull; ++minK) {
            INFO("minK " << minK);
            for (size_t maxK{minK}; maxK < 5ull; ++maxK) {
                INFO("maxK " << maxK);
                CHECK(ss::isComplete(gen::bestKnown(N, minK, maxK), minK, maxK));
            }
        }
    }
}

TEST_CASE("check search scheme generator h2 for completness", "[isComplete][h2]") {
    // Note N must be larger than maxK
    for (size_t N{1}; N < 10ull; ++N) { // Number of pieces
        INFO("N " << N);
        for (size_t minK{0}; minK < std::min(N, size_t{5}); ++minK) {
            INFO("minK " << minK);
            for (size_t maxK{minK}; maxK < std::min(N, size_t{5}); ++maxK) {
                INFO("maxK " << maxK);
                CHECK(ss::isComplete(gen::h2(N, minK, maxK), minK, maxK));
            }
        }
    }
}

TEST_CASE("check search scheme generator kianfar for completness", "[isComplete][kianfar]") {
    // Kianfar only exists for certain k's
    for (size_t k{0}; k < 4; ++k) {
        CHECK(ss::isComplete(gen::kianfar(k), 0, k));
    }
}

TEST_CASE("check search scheme generator kucherov for completness", "[isComplete][kucherov]") {
    // Kucherov only exists for certain N and k combinations
    for (size_t k{0}; k < 5; ++k) {
        CHECK(ss::isComplete(gen::kucherov(k+1, k), 0, k));
    }
    for (size_t k{0}; k < 5; ++k) {
        CHECK(ss::isComplete(gen::kucherov(k+2, k), 0, k));
    }
}

TEST_CASE("check search scheme generator optimum for completness", "[isComplete][optimum]") {
    for (size_t minK{0}; minK < 4; ++minK) {
        for (size_t maxK{minK}; maxK < 4; ++maxK) {
            CHECK(ss::isComplete(gen::optimum(minK, maxK), minK, maxK));
        }
    }
}

TEST_CASE("check search scheme generator pigeon for completness", "[isComplete][pigeon]") {
    for (size_t minK{0}; minK < 10; ++minK) {
        for (size_t maxK{minK}; maxK < 10; ++maxK) {
            CHECK(ss::isComplete(gen::pigeon_trivial(minK, maxK), minK, maxK));
        }
    }
    for (size_t minK{0}; minK < 10; ++minK) {
        for (size_t maxK{minK}; maxK < 10; ++maxK) {
            CHECK(ss::isComplete(gen::pigeon_opt(minK, maxK), minK, maxK));
        }
    }
}

TEST_CASE("check search scheme generator suffixFilter for completness", "[isComplete][suffixFilter]") {
    // Note N must be larger than maxK
    for (size_t N{1}; N < 10; ++N) { // Number of pieces
        INFO("N " << N);
        for (size_t minK{0}; minK < std::min(N, size_t{5}); ++minK) {
            INFO("minK " << minK);
            for (size_t maxK{minK}; maxK < std::min(N, size_t{5}); ++maxK) {
                INFO("maxK " << maxK);
                CHECK(ss::isComplete(gen::suffixFilter(N, minK, maxK), minK, maxK));
            }
        }
    }
}

TEST_CASE("check search scheme generator 01*0 (zeroOnesZero) for completness", "[isComplete][zeroOnesZero]") {
    for (size_t minK{0}; minK < 10; ++minK) {
        for (size_t maxK{minK}; maxK < 10; ++maxK) {
            CHECK(ss::isComplete(gen::zeroOnesZero_trivial(minK, maxK), minK, maxK));
        }
    }
    for (size_t minK{0}; minK < 10; ++minK) {
        for (size_t maxK{minK}; maxK < 10; ++maxK) {
            CHECK(ss::isComplete(gen::zeroOnesZero_opt(minK, maxK), minK, maxK));
        }
    }
}
