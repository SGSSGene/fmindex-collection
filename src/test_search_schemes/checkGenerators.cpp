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

TEST_CASE("check search scheme generator backtracking", "[isValid][backtracking]") {
    for (size_t N{1}; N < 20; ++N) { // Number of pieces
        INFO("N " << N);
        for (size_t minK{0}; minK < 10ull; ++minK) {
            INFO("minK " << minK);
            for (size_t maxK{minK}; maxK < 10ull; ++maxK) {
                INFO("maxK " << maxK);
                CHECK(ss::isValid(gen::backtracking(N, minK, maxK)));
            }
        }
    }
}

TEST_CASE("check search scheme generator bestKnown", "[isValid][bestKnown]") {
    for (size_t N{1}; N < 20; ++N) { // Number of pieces
        INFO("N " << N);
        for (size_t minK{0}; minK < 10ull; ++minK) {
            INFO("minK " << minK);
            for (size_t maxK{minK}; maxK < 10ull; ++maxK) {
                INFO("maxK " << maxK);
                CHECK(ss::isValid(gen::bestKnown(N, minK, maxK)));
            }
        }
    }
}

TEST_CASE("check search scheme generator h2", "[isValid][h2]") {
    // Note N must be larger than maxK
    for (size_t N{1}; N < 20; ++N) { // Number of pieces
        INFO("N " << N);
        for (size_t minK{0}; minK < std::min(N, size_t{10}); ++minK) {
            INFO("minK " << minK);
            for (size_t maxK{minK}; maxK < std::min(N, size_t{10}); ++maxK) {
                INFO("maxK " << maxK);
                CHECK(ss::isValid(gen::h2(N, minK, maxK)));
            }
        }
    }
}

TEST_CASE("check search scheme generator kianfar", "[isValid][kianfar]") {
    // Kianfar only exists for certain k's
    for (size_t k{0}; k < 4; ++k) {
        CHECK(ss::isValid(gen::kianfar(k)));
    }
}

TEST_CASE("check search scheme generator kucherov", "[isValid][kucherov]") {
    // Kucherov only exists for certain N and k combinations
    for (size_t k{0}; k < 5; ++k) {
        CHECK(ss::isValid(gen::kucherov(k+1, k)));
    }
    for (size_t k{0}; k < 5; ++k) {
        CHECK(ss::isValid(gen::kucherov(k+2, k)));
    }
}

TEST_CASE("check search scheme generator optimum", "[isValid][optimum]") {
    for (size_t minK{0}; minK < 4; ++minK) {
        for (size_t maxK{minK}; maxK < 4; ++maxK) {
            CHECK(ss::isValid(gen::optimum(minK, maxK)));
        }
    }
}

TEST_CASE("check search scheme generator pigeon", "[isValid][pigeon]") {
    for (size_t minK{0}; minK < 20; ++minK) {
        for (size_t maxK{minK}; maxK < 20; ++maxK) {
            CHECK(ss::isValid(gen::pigeon_trivial(minK, maxK)));
        }
    }
    for (size_t minK{0}; minK < 20; ++minK) {
        for (size_t maxK{minK}; maxK < 20; ++maxK) {
            CHECK(ss::isValid(gen::pigeon_opt(minK, maxK)));
        }
    }
}

TEST_CASE("check search scheme generator suffixFilter", "[isValid][suffixFilter]") {
    // Note N must be larger than maxK
    for (size_t N{1}; N < 20; ++N) { // Number of pieces
        INFO("N " << N);
        for (size_t minK{0}; minK < std::min(N, size_t{10}); ++minK) {
            INFO("minK " << minK);
            for (size_t maxK{minK}; maxK < std::min(N, size_t{10}); ++maxK) {
                INFO("maxK " << maxK);
                CHECK(ss::isValid(gen::suffixFilter(N, minK, maxK)));
            }
        }
    }
}

TEST_CASE("check search scheme generator 01*0 (zeroOnesZero)", "[isValid][zeroOnesZero]") {
    for (size_t minK{0}; minK < 20; ++minK) {
        for (size_t maxK{minK}; maxK < 20; ++maxK) {
            CHECK(ss::isValid(gen::zeroOnesZero_trivial(minK, maxK)));
        }
    }
    for (size_t minK{0}; minK < 20; ++minK) {
        for (size_t maxK{minK}; maxK < 20; ++maxK) {
            CHECK(ss::isValid(gen::zeroOnesZero_opt(minK, maxK)));
        }
    }
}
