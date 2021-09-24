#include <search_schemes/isValid.h>
#include <catch2/catch.hpp>

namespace ss = search_schemes;

TEST_CASE("check is valid", "[isValid]") {

    auto search = ss::Search{
        {0},
        {0},
        {0}
    };
    REQUIRE(ss::isValid(search));

}
