#include <search_schemes/isComplete.h>
#include <catch2/catch.hpp>

namespace ss = search_schemes;

TEST_CASE("check is complete", "[isComplete]") {

    CHECK(ss::isComplete(ss::Scheme{ss::Search {
        {0, 1},
        {0, 0},
        {0, 0},
    }}, 0, 0));

    CHECK(not ss::isComplete(ss::Scheme{ss::Search {
        {0, 1},
        {0, 0},
        {0, 0},
    }}, 0, 1));

    CHECK(ss::isComplete(ss::Scheme{ss::Search {
        {0, 1},
        {0, 0},
        {1, 1},
    }}, 0, 1));

    CHECK(ss::isComplete(ss::Scheme{ss::Search {
        {0, 1},
        {0, 1},
        {1, 1},
    }}, 1, 1));
}
