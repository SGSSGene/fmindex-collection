#include <search_schemes/nodeCount.h>
#include <search_schemes/expand.h>
#include <search_schemes/generator/backtracking.h>
#include <catch2/catch.hpp>

namespace ss = search_schemes;
namespace gen = ss::generator;

TEST_CASE("check node counts", "[nodeCount]") {
    SECTION("known length of search schemes with 0 errors") {
        for (size_t n{1}; n < 1000; ++n) {
            INFO("length(n): " << n);
            CHECK(n == ss::nodeCount(gen::backtracking(n, 0, 0), 4));
        }

        for (size_t n{1}; n < 1000; ++n) {
            INFO("length(n): " << n);
            CHECK(n == ss::nodeCount(ss::expand(gen::backtracking(1, 0, 0), n), 4));
        }
    }

    SECTION("known length of some search schemes") {
        CHECK( 4 == ss::nodeCount(gen::backtracking(1, 0, 1), 4));
        CHECK(11 == ss::nodeCount(gen::backtracking(2, 0, 1), 4));
        CHECK(21 == ss::nodeCount(gen::backtracking(3, 0, 1), 4));
        CHECK(20 == ss::nodeCount(gen::backtracking(2, 0, 2), 4));
    }


}

