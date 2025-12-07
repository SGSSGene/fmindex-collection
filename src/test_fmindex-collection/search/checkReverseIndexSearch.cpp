// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include "../string/allStrings.h"

#include <catch2/catch_all.hpp>
#include <fmindex-collection/fmindex/ReverseFMIndex.h>
#include <fmindex-collection/search/BacktrackingWithBuffers.h>
#include <fmindex-collection/search_scheme/generator/all.h>

TEST_CASE("searching with collection and backtracking with buffers on a reversed fmindex", "[collection][search][reverse]") {
    using String = fmc::string::InterleavedBitvector16<255>;

    auto input  = std::vector<std::vector<uint8_t>>{{'A', 'A', 'A', 'B', 'A', 'A', 'A', 'C', 'A', 'A', 'A'},
                                                    {'A', 'A', 'A', 'C', 'A', 'A', 'A', 'B', 'A', 'A', 'A'}};

    auto index = fmc::ReverseFMIndex<String>{input, /*samplingRate*/1, /*threadNbr*/1};

    auto query = std::vector<uint8_t>{'A', 'C'};
    auto search_scheme = fmc::search_scheme::generator::backtracking(1, 0, 0);

    using cursor_t = fmc::select_cursor_t<decltype(index)>;

    auto buffer1 = std::vector<std::pair<cursor_t, size_t>>{};
    auto buffer2 = std::vector<std::pair<cursor_t, size_t>>{};

    fmc::search_backtracking_with_buffers::search(index, query, /*error*/ 0, buffer1, buffer2, [&](auto cursor, auto errors) {
        CHECK(errors == 0);
        CHECK(cursor.lb == 22);
        CHECK(cursor.count() == 2);

        auto result = std::vector<std::pair<size_t, size_t>>{};
        for (size_t i{0}; i < cursor.count(); ++i) {
            auto [il, pl, offset] = index.locate(i + cursor.lb);
            pl += offset;
            result.emplace_back(il, pl);
        }
        auto expected = std::vector<std::pair<size_t, size_t>>{{1, 4}, {0, 8}};
        REQUIRE(result.size() == expected.size());
        for (size_t i{0}; i < result.size(); ++i) {
            CHECK(std::get<0>(result[i]) == std::get<0>(expected[i]));
            CHECK(std::get<1>(result[i]) == std::get<1>(expected[i]));
        }
    });
}
