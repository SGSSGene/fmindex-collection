// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <catch2/catch_all.hpp>
#include <fmindex-collection/search_scheme/generator/backtracking.h>
#include <fmindex-collection/search_scheme/generator/bestKnown.h>
#include <fmindex-collection/search_scheme/generator/h2.h>
#include <fmindex-collection/search_scheme/generator/hato.h>
#include <fmindex-collection/search_scheme/generator/kianfar.h>
#include <fmindex-collection/search_scheme/generator/kucherov.h>
#include <fmindex-collection/search_scheme/generator/optimum.h>
#include <fmindex-collection/search_scheme/generator/pigeon.h>
#include <fmindex-collection/search_scheme/generator/suffixFilter.h>
#include <fmindex-collection/search_scheme/generator/zeroOnesZero.h>
#include <fmindex-collection/search_scheme/isComplete.h>
#include <fmindex-collection/search_scheme/isValid.h>

namespace ss = fmc::search_scheme;
namespace gen = ss::generator;

TEST_CASE("check search scheme generator backtracking for completeness", "[isComplete][backtracking]") {
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

TEST_CASE("check search scheme generator bestKnown for completeness", "[isComplete][bestKnown]") {
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

TEST_CASE("check search scheme generator h2 for completeness", "[isComplete][h2]") {
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

TEST_CASE("check search scheme generator hato for completeness", "[isComplete][hato]") {
    for (size_t maxK{0}; maxK < 8; ++maxK) {
        INFO("maxK " << maxK);
        CHECK(ss::isComplete(gen::hato(maxK), 0, maxK));
    }
}

TEST_CASE("check search scheme generator kianfar for completeness", "[isComplete][kianfar]") {
    // Kianfar only exists for certain k's
    for (size_t k{0}; k < 4; ++k) {
        CHECK(ss::isComplete(gen::kianfar(k), 0, k));
    }
}

TEST_CASE("check search scheme generator kucherov for completeness", "[isComplete][kucherov]") {
    // Kucherov only exists for certain N and k combinations
    for (size_t k{0}; k < 5; ++k) {
        CHECK(ss::isComplete(gen::kucherov(k+1, k), 0, k));
    }
    for (size_t k{0}; k < 5; ++k) {
        CHECK(ss::isComplete(gen::kucherov(k+2, k), 0, k));
    }
}

TEST_CASE("check search scheme generator optimum for completeness", "[isComplete][optimum]") {
    for (size_t minK{0}; minK < 4; ++minK) {
        for (size_t maxK{minK}; maxK < 4; ++maxK) {
            CHECK(ss::isComplete(gen::optimum(minK, maxK), minK, maxK));
        }
    }
}

TEST_CASE("check search scheme generator pigeon for completeness", "[isComplete][pigeon]") {
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

TEST_CASE("check search scheme generator suffixFilter for completeness", "[isComplete][suffixFilter]") {
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

TEST_CASE("check search scheme generator 01*0 (zeroOnesZero) for completeness", "[isComplete][zeroOnesZero]") {
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
