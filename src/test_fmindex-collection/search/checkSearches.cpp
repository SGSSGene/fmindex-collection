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


TEST_CASE("check searches with errors", "[searches]") {
    using OccTable = fmindex_collection::occtable::EprV2_16<256>;
    using Index = fmindex_collection::BiFMIndex<OccTable>;

    auto input  = std::vector<std::vector<uint8_t>>{{'A', 'A', 'A', 'C', 'A', 'A', 'A', 'B', 'A', 'A', 'A'},
                                                    {'A', 'A', 'A', 'B', 'A', 'A', 'A', 'C', 'A', 'A', 'A'}};

    auto index = Index{input, /*samplingRate*/1, /*threadNbr*/1};

    auto queries = std::vector<std::vector<uint8_t>> {std::vector<uint8_t>{'C', 'C'}, std::vector<uint8_t>{'B', 'B'}};

    SECTION("backtracking, single search") {
        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        for (size_t qidx{0}; qidx != queries.size(); ++qidx) {
            auto const& query = queries[qidx];
            fmindex_collection::search_backtracking::search(index, query, 1, [&](auto cursor, auto errors) {
                (void)errors;
                for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                    results.emplace_back(qidx, sid, spos);
                }
            });
        }

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 2},
            {0, 0, 3},
            {0, 1, 6},
            {0, 1, 7},
            {1, 0, 6},
            {1, 0, 7},
            {1, 1, 2},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("backtracking, all search") {
        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_backtracking::search(index, queries, 1, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 2},
            {0, 0, 3},
            {0, 1, 6},
            {0, 1, 7},
            {1, 0, 6},
            {1, 0, 7},
            {1, 1, 2},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("backtracking with buffers") {
        auto buffer1 = std::vector<std::pair<fmindex_collection::select_cursor_t<Index>, size_t>>{};
        auto buffer2 = std::vector<std::pair<fmindex_collection::select_cursor_t<Index>, size_t>>{};

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        for (size_t qidx{0}; qidx < queries.size(); ++qidx) {
            auto const& query = queries[qidx];
            fmindex_collection::search_backtracking_with_buffers::search(index, query, 1, buffer1, buffer2, [&](auto cursor, auto errors) {
                (void)errors;
                for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                    results.emplace_back(qidx, sid, spos);
                }
            });
        }

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 2},
            {0, 0, 3},
            {0, 1, 6},
            {0, 1, 7},
            {1, 0, 6},
            {1, 0, 7},
            {1, 1, 2},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("search no errors, all search") {
        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_no_errors::search(index, queries, [&](auto qidx, auto cursor) {
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
        };
        CHECK(results == expected);
    }

    SECTION("search one error, all search") {
        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_one_error::search(index, queries, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 2},
            {0, 0, 3},
            {0, 1, 6},
            {0, 1, 7},
            {1, 0, 6},
            {1, 0, 7},
            {1, 1, 2},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("pseudo search, all search") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_pseudo::search</*EditDistance=*/true>(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
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

    SECTION("pseudo search, single searches") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        for (size_t qidx{0}; qidx < queries.size(); ++qidx) {
            auto const& query = queries[qidx];
            fmindex_collection::search_pseudo::search</*EditDistance=*/true>(index, query, search_scheme, [&](auto cursor, auto errors) {
                (void)errors;
                for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                    results.emplace_back(qidx, sid, spos);
                }
            });
        }

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


    SECTION("pseudo search, hamming distance, all search") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_pseudo::search</*EditDistance=*/false>(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 2},
            {0, 0, 3},
            {0, 1, 6},
            {0, 1, 7},
            {1, 0, 6},
            {1, 0, 7},
            {1, 1, 2},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }


    SECTION("search ng12, all search") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_ng12::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
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

    SECTION("search ng14, all search") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_ng14::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 0, 3},
            {0, 1, 7},
            {0, 1, 7},
            {1, 0, 7},
            {1, 0, 7},
            {1, 1, 3},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("search ng15, all search") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_ng15::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 2},
            {0, 0, 3},
            {0, 1, 6},
            {0, 1, 7},
            {1, 0, 6},
            {1, 0, 7},
            {1, 1, 2},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }
    //!TODO Doesn't work for such short search schemes??
    SECTION("search ng16, all search") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::backtracking(2, 0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_ng16::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
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

    SECTION("search ng17, all search") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_ng17::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 2},
            {0, 0, 3},
            {0, 1, 6},
            {0, 1, 7},
            {1, 0, 6},
            {1, 0, 7},
            {1, 1, 2},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

/*    SECTION("search ng20, all search") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::backtracking(2, 0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_ng20::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
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
    }*/


    SECTION("search ng21, all search") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_ng21::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 0, 3},
            {0, 1, 7},
            {0, 1, 7},
            {1, 0, 7},
            {1, 0, 7},
            {1, 1, 3},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("search ng21, all search_n") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_ng21::search_n(index, queries, search_scheme, 3, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 1, 7},
            {0, 1, 7},
            {1, 0, 7},
            {1, 0, 7},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("search ng21, all search_best") {
        auto search_scheme_e0 = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 0), queries[0].size());
        auto search_scheme_e1 = search_schemes::expand(search_schemes::generator::pigeon_opt(1, 1), queries[0].size());
        auto search_scheme_e2 = search_schemes::expand(search_schemes::generator::pigeon_opt(2, 2), queries[0].size());
        auto search_schemes = std::vector{search_scheme_e0, search_scheme_e1, search_scheme_e2};

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_ng21::search_best(index, queries, search_schemes, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 0, 3},
            {0, 1, 7},
            {0, 1, 7},
            {1, 0, 7},
            {1, 0, 7},
            {1, 1, 3},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("search ng21, all search_best_n") {
        auto search_scheme_e0 = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 0), queries[0].size());
        auto search_scheme_e1 = search_schemes::expand(search_schemes::generator::pigeon_opt(1, 1), queries[0].size());
        auto search_schemes = std::vector{search_scheme_e0, search_scheme_e1};

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_ng21::search_best_n(index, queries, search_schemes, 3, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 1, 7},
            {0, 1, 7},
            {1, 0, 7},
            {1, 0, 7},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("search ng21 V2, all search") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_ng21V2::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 0, 3},
            {0, 1, 7},
            {0, 1, 7},
            {1, 0, 7},
            {1, 0, 7},
            {1, 1, 3},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("search ng21 V3, all search") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_ng21V3::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 0, 3},
            {0, 1, 7},
            {0, 1, 7},
            {1, 0, 7},
            {1, 0, 7},
            {1, 1, 3},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("search ng21 V4, all search") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_ng21V4::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 0, 3},
            {0, 1, 7},
            {0, 1, 7},
            {1, 0, 7},
            {1, 0, 7},
            {1, 1, 3},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("search ng21 V5, all search") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_ng21V5::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 0, 3},
            {0, 1, 7},
            {0, 1, 7},
            {1, 0, 7},
            {1, 0, 7},
            {1, 1, 3},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("search ng21 V6, all search") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_ng21V6::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 0, 3},
            {0, 1, 7},
            {0, 1, 7},
            {1, 0, 7},
            {1, 0, 7},
            {1, 1, 3},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("search ng21 V6, all search_n") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_ng21V6::search_n(index, queries, search_scheme, 3, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 1, 7},
            {0, 1, 7},
            {1, 0, 7},
            {1, 0, 7},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("search ng21 V6, all search_best") {
        auto search_scheme_e0 = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 0), queries[0].size());
        auto search_scheme_e1 = search_schemes::expand(search_schemes::generator::pigeon_opt(1, 1), queries[0].size());
        auto search_scheme_e2 = search_schemes::expand(search_schemes::generator::pigeon_opt(2, 2), queries[0].size());
        auto search_schemes = std::vector{search_scheme_e0, search_scheme_e1, search_scheme_e2};

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_ng21V6::search_best(index, queries, search_schemes, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 0, 3},
            {0, 1, 7},
            {0, 1, 7},
            {1, 0, 7},
            {1, 0, 7},
            {1, 1, 3},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("search ng21 V6, all search_best_n") {
        auto search_scheme_e0 = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 0), queries[0].size());
        auto search_scheme_e1 = search_schemes::expand(search_schemes::generator::pigeon_opt(1, 1), queries[0].size());
        auto search_schemes = std::vector{search_scheme_e0, search_scheme_e1};

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_ng21V6::search_best_n(index, queries, search_schemes, 3, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 1, 7},
            {0, 1, 7},
            {1, 0, 7},
            {1, 0, 7},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("search ng21 V7, all search") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_ng21V7::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 0, 3},
            {0, 1, 7},
            {0, 1, 7},
            {1, 0, 7},
            {1, 0, 7},
            {1, 1, 3},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("search ng21 V7, all search_n") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_ng21V7::search_n(index, queries, search_scheme, 3, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 1, 7},
            {0, 1, 7},
            {1, 0, 7},
            {1, 0, 7},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("search ng21 V7, all search_best") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 2), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_ng21V7::search_best(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 0, 3},
            {0, 1, 7},
            {0, 1, 7},
            {1, 0, 7},
            {1, 0, 7},
            {1, 1, 3},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("search ng21 V7, all search_best_n") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 2), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_ng21V7::search_best_n(index, queries, search_scheme, 3, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 1, 7},
            {0, 1, 7},
            {1, 0, 7},
            {1, 0, 7},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }


    SECTION("search ng22, all search") {
        auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmindex_collection::search_ng22::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors, auto const& action) {
            (void)errors;
            (void)action;
            for (auto [sid, spos] : fmindex_collection::LocateLinear{index, cursor}) {
                results.emplace_back(qidx, sid, spos);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 2},
            {0, 0, 3},
            {0, 1, 6},
            {0, 1, 7},
            {1, 0, 6},
            {1, 0, 7},
            {1, 1, 2},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }
}
