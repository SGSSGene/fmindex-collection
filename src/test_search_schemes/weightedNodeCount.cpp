// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include <catch2/catch_all.hpp>
#include <search_schemes/expand.h>
#include <search_schemes/generator/backtracking.h>
#include <search_schemes/weightedNodeCount.h>

namespace ss = search_schemes;
namespace gen = ss::generator;

TEST_CASE("check node counts", "[weightedNodeCount]") {
    SECTION("known length of search schemes with 0 errors") {
        for (size_t n{1}; n < 14; ++n) {
            INFO("length(n): " << n);
            CHECK(n == ss::weightedNodeCount</*Edit=*/false>(gen::backtracking(n, 0, 0), 4, 1'000'000'000));
        }

        for (size_t n{1}; n < 14; ++n) {
            INFO("length(n): " << n);
            CHECK(n == ss::weightedNodeCount</*Edit=*/false>(ss::expand(gen::backtracking(1, 0, 0), n), 4, 1'000'000'000));
        }

        for (size_t n{15}; n < 1000; ++n) {
            INFO("length(n): " << n);
            CHECK(16 > ss::weightedNodeCount</*Edit=*/false>(gen::backtracking(n, 0, 0), 4, 1'000'000'000));
        }

        for (size_t n{15}; n < 1000; ++n) {
            INFO("length(n): " << n);
            CHECK(16 >  ss::weightedNodeCount</*Edit=*/false>(ss::expand(gen::backtracking(1, 0, 0), n), 4, 1'000'000'000));
        }

    }

    SECTION("known length of some search schemes") {
        CHECK( 4 == ss::weightedNodeCount</*Edit=*/false>(gen::backtracking(1, 0, 1), 4, 1'000'000'000));
        CHECK(11 == ss::weightedNodeCount</*Edit=*/false>(gen::backtracking(2, 0, 1), 4, 1'000'000'000));
        CHECK(21 == ss::weightedNodeCount</*Edit=*/false>(gen::backtracking(3, 0, 1), 4, 1'000'000'000));
        CHECK(20 == ss::weightedNodeCount</*Edit=*/false>(gen::backtracking(2, 0, 2), 4, 1'000'000'000));
    }


}

