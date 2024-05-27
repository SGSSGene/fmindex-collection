// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <catch2/catch_all.hpp>
#include <fmindex-collection/fmindex/BiFMIndex.h>
#include <fmindex-collection/locate.h>
#include <fmindex-collection/occtable/all.h>
#include <fmindex-collection/search/all.h>
#include <search_schemes/generator/all.h>
#include <search_schemes/expand.h>


TEST_CASE("locating using LocateFMTree", "[locate][fmtree]") {
    using OccTable = fmindex_collection::occtable::EprV2_16<256>;
    using Index = fmindex_collection::BiFMIndex<OccTable>;

    auto input  = std::vector<std::vector<uint8_t>>{{'A', 'A', 'A', 'C', 'A', 'A', 'A', 'B', 'A', 'A', 'A'},
                                                    {'A', 'A', 'A', 'B', 'A', 'A', 'A', 'C', 'A', 'A', 'A'}};

    size_t samplingRate = 4;
    auto index = Index{input, samplingRate, /*threadNbr*/1};

    auto queries = std::vector<std::vector<uint8_t>> {std::vector<uint8_t>{'C', 'C'}, std::vector<uint8_t>{'B', 'B'}};

    SECTION("test LocateFMTree") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_pseudo::search</*EditDistance=*/true>(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateFMTree{index, cursor, samplingRate, /*.maxDepth=*/3}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 2},
            {0, 0, 3},
            {0, 0, 3},
            {0, 0, 3},
            {0, 1, 6},
            {0, 1, 7},
            {0, 1, 7},
            {0, 1, 7},
            {1, 0, 6},
            {1, 0, 7},
            {1, 0, 7},
            {1, 0, 7},
            {1, 1, 2},
            {1, 1, 3},
            {1, 1, 3},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("test locateFMTree") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_pseudo::search</*EditDistance=*/true>(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            fmindex_collection::locateFMTree<4>(index, cursor, [&](size_t sid, size_t spos) {
                results.emplace_back(qidx, sid, spos);
            }, samplingRate);
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 2},
            {0, 0, 3},
            {0, 0, 3},
            {0, 0, 3},
            {0, 1, 6},
            {0, 1, 7},
            {0, 1, 7},
            {0, 1, 7},
            {1, 0, 6},
            {1, 0, 7},
            {1, 0, 7},
            {1, 0, 7},
            {1, 1, 2},
            {1, 1, 3},
            {1, 1, 3},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

}
