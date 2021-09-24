#include <search_schemes/isValid.h>
#include <search_schemes/generator/backtracking.h>
#include <search_schemes/generator/bestKnown.h>
#include <search_schemes/generator/h2.h>
#include <search_schemes/generator/kianfar.h>
#include <search_schemes/generator/kucherov.h>
#include <search_schemes/generator/optimum.h>
#include <search_schemes/generator/pigeon.h>
#include <catch2/catch.hpp>

namespace ss = search_schemes;
namespace gen = ss::generator;

TEST_CASE("check search scheme generator backtracking", "[isValid][backtracking]") {
    for (size_t N{1}; N < 10; ++N) { // Number of pieces
        INFO("N " << N);
        for (size_t minK{0}; minK < 6; ++minK) {
            INFO("minK " << minK);
            for (size_t maxK{minK}; maxK < 6; ++maxK) {
                INFO("maxK " << maxK);
                CHECK(ss::isValid(gen::backtracking(N, minK, maxK)));
            }
        }
    }
}

TEST_CASE("check search scheme generator bestKnown", "[isValid][bestKnown]") {
    for (size_t N{1}; N < 10; ++N) { // Number of pieces
        INFO("N " << N);
        for (size_t minK{0}; minK < 6; ++minK) {
            INFO("minK " << minK);
            for (size_t maxK{minK}; maxK < 6; ++maxK) {
                INFO("maxK " << maxK);
                CHECK(ss::isValid(gen::bestKnown(N, minK, maxK)));
            }
        }
    }
}

TEST_CASE("check search scheme generator h2", "[isValid][h2]") {
    // Note N must be larger than maxK
    for (size_t N{1}; N < 10; ++N) { // Number of pieces
        INFO("N " << N);
        for (size_t minK{0}; minK < std::min(N, 6ul); ++minK) {
            INFO("minK " << minK);
            for (size_t maxK{minK}; maxK < std::min(N, 6ul); ++maxK) {
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
