#include <search_schemes/isValid.h>
#include <search_schemes/generator/backtracking.h>
#include <search_schemes/generator/bestKnown.h>
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
