// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <catch2/catch_all.hpp>
#include <fmindex-collection/fmindex/BiFMIndex.h>
#include <fmindex-collection/locate.h>
#include <fmindex-collection/occtable/all.h>
#include <fmindex-collection/search/SearchHammingSM.h>
#include <fmindex-collection/search/all.h>
#include <fmindex-collection/search_scheme/generator/all.h>
#include <fmindex-collection/search_scheme/expand.h>
#include <nanobench.h>


TEST_CASE("check search with hamming and scoring matrix with errors", "[searches][hamming][scoring-matrix]") {
//    using OccTable = fmindex_collection::occtable::EprV7<21>;
    using OccTable = fmindex_collection::occtable::Interleaved_16<21>;
    using Index = fmindex_collection::BiFMIndex<OccTable>;

    SECTION("hamming_sm, all search") {
        srand(0);
        auto input = std::vector<std::vector<uint8_t>>{{20}};
        auto& ref = input.back();

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>>{};
        for (size_t i = 0; i < 10000000; ++i) {
        //for (size_t i = 0; i < 100; ++i) {
            ref.emplace_back(rand()%20+1);
        }
        ref.emplace_back(20);

        auto queries = std::vector<std::vector<uint8_t>> {};
        constexpr size_t QueryLen = 16;
        for (size_t i{0}; i < 1000; ++i) {
        //for (size_t i{0}; i < 100; ++i) {
            auto p = input[0].size()-QueryLen;
            auto pos = (rand()+1) % p;
            queries.emplace_back();
            for (size_t j{}; j < QueryLen; ++j) {
                queries.back().push_back(input[0][pos+j]);
            }
            expected.emplace_back(queries.size()-1, 0, pos);
            // introduce an error
            {
                auto pos = rand() % QueryLen;
                auto inc = rand() % 19 + 1;
                queries.back()[pos] = ((queries.back()[pos] + inc) % 20) + 1;
            }
        }
/*        expected.erase(expected.begin(), expected.begin()+97);
        queries.erase(queries.begin(), queries.begin()+97);
        expected.erase(expected.end()-2, expected.end());
        queries.erase(queries.end()-2, queries.end());
        for (size_t i{0}; i < expected.size(); ++i) {
            auto& [qid, sid, spos] = expected[i];
            qid = i;
        }*/

/*        std::cout << "refs:\n";
        for (auto r : input) {
            for (auto v : r) {
                std::cout << (char)(v+'A');
            }
            std::cout << "\n";
        }
        std::cout << "\n\n";
        std::cout << "queries:\n";
        for (auto r : queries) {
            for (auto v : r) {
                std::cout << (char)(v+'A');
            }
            std::cout << "\n";
        }
        std::cout << "\n\n";*/

        auto index = Index{input, /*samplingRate*/1, /*threadNbr*/1};


        auto search_scheme = fmindex_collection::search_scheme::generator::pigeon_opt(0, 1);

        {
            auto sm = fmindex_collection::search_hamming_sm::ScoringMatrix<28, 21>{}; // Inversed-Identity is set by default
            sm.setCost(21,  5, 0);
            sm.setCost(22, 13, 0);
            sm.setCost(23,  4, 0);
            sm.setCost(24,  7, 0);
            sm.setCost(25, 12, 0);
            sm.setCost(26, 17, 0);
            sm.setCost(27, 19, 0);

            ankerl::nanobench::Bench().minEpochTime(std::chrono::milliseconds{100}).run("hamming with special scoring matrix and reducing alphabet", [&]() {
                fmindex_collection::search_hamming_sm::search(index, queries, search_scheme, sm, [&](auto qidx, auto cursor, auto errors) {
                    (void)errors;
                    (void)qidx;
                    ankerl::nanobench::doNotOptimizeAway(cursor);
                });
            });

            auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
            fmindex_collection::search_hamming_sm::search(index, queries, search_scheme, sm, [&](auto qidx, auto cursor, auto errors) {
                (void)errors;
                (void)qidx;
                for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                    results.emplace_back(qidx, sid, spos);
                }
            });
            std::ranges::sort(results);
            std::ranges::sort(expected);

            CHECK(results.size() == expected.size());
            CHECK(results == expected);
        }

        {
            auto expanded_search_scheme = fmindex_collection::search_scheme::expand(search_scheme, queries[0].size());

            ankerl::nanobench::Bench().minEpochTime(std::chrono::milliseconds{10}).run("hamming without scoring matrix", [&]() {
                fmindex_collection::search_pseudo::search</*EditDistance=*/false>(index, queries, expanded_search_scheme, [&](auto qidx, auto cursor, auto errors) {
                    (void)errors;
                    (void)qidx;
                    ankerl::nanobench::doNotOptimizeAway(cursor);
                });
            });

            auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
            fmindex_collection::search_pseudo::search</*EditDistance=*/false>(index, queries, expanded_search_scheme, [&](auto qidx, auto cursor, auto errors) {
                (void)errors;
                (void)qidx;
                for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                    results.emplace_back(qidx, sid, spos);
                }
            });
            std::ranges::sort(results);
            std::ranges::sort(expected);

            CHECK(results == expected);
        }
    }
}
