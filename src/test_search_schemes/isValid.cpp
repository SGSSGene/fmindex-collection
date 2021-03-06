#include <search_schemes/isValid.h>
#include <catch2/catch.hpp>

namespace ss = search_schemes;

TEST_CASE("check is valid", "[isValid]") {

    auto search = ss::Search{
        {0},
        {0},
        {0}
    };
    CHECK(ss::isValid(search));

    CHECK(ss::isValid(ss::Search {
        {0, 1},
        {0, 0},
        {0, 0},
    }));

    CHECK(ss::isValid(ss::Search {
        {1, 0},
        {0, 0},
        {0, 0},
    }));
    CHECK(ss::isValid(ss::Search {
        {0, 1, 2},
        {0, 0, 0},
        {0, 0, 0},
    }));
    CHECK(ss::isValid(ss::Search {
        {1, 0, 2},
        {0, 0, 0},
        {0, 0, 0},
    }));
    CHECK(ss::isValid(ss::Search {
        {1, 2, 0},
        {0, 0, 0},
        {0, 0, 0},
    }));
    CHECK(ss::isValid(ss::Search {
        {2, 1, 0},
        {0, 0, 0},
        {0, 0, 0},
    }));
    CHECK(not ss::isValid(ss::Search {
        {0, 2, 1},
        {0, 0, 0},
        {0, 0, 0},
    }));
    CHECK(not ss::isValid(ss::Search {
        {2, 0, 1},
        {0, 0, 0},
        {0, 0, 0},
    }));
    CHECK(not ss::isValid(ss::Search {
        {0, 0, 2},
        {0, 0, 0},
        {0, 0, 0},
    }));
}
