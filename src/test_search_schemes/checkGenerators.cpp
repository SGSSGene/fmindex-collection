#include <search_schemes/isValid.h>
#include <search_schemes/generator/backtracking.h>
#include <catch2/catch.hpp>

namespace ss = search_schemes;
namespace gen = ss::generator;

TEST_CASE("check search scheme generator backtracking", "[isValid][backtracking]") {
    for (size_t N{1}; N < 10; ++N) { // Number of pieces
        INFO("N " << N);
        for (size_t minK{0}; minK < 3; ++minK) {
            INFO("minK " << minK);
            for (size_t maxK{minK}; maxK < 4; ++maxK) {
                INFO("maxK " << maxK);
                CHECK(ss::isValid(gen::backtracking(N, minK, maxK)));
            }
        }
    }

}
