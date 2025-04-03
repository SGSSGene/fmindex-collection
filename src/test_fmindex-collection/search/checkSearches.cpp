// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <catch2/catch_all.hpp>
#include <fmindex-collection/fmindex/BiFMIndex.h>
#include <fmindex-collection/locate.h>
#include <fmindex-collection/occtable/all.h>
#include <fmindex-collection/search/all.h>
#include <nanobench.h>
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

TEST_CASE("benchmark searches with errors", "[searches][!benchmark]") {
    SECTION("benchmarking") {
        using OccTable = fmindex_collection::occtable::EprV2_16<256>;
        using Index = fmindex_collection::BiFMIndex<OccTable>;

        srand(0);

        auto generateSequence = [](size_t l) {
            auto seq = std::vector<uint8_t>{};
            seq.reserve(l);
            for (size_t i{0}; i < l; ++i) {
                seq.emplace_back(1 + rand()%4);
            }
            return seq;
        };

        auto generateSequences = [generateSequence](size_t n, size_t l) {
            auto ref = std::vector<std::vector<uint8_t>>{};
            ref.reserve(n);
            for (size_t i{0}; i < n; ++i) {
                ref.emplace_back(generateSequence(l));
            }
            return ref;
        };

        auto generateRead = [](std::vector<std::vector<uint8_t>> const& ref, size_t l, size_t e) {
            auto subjId = rand() % ref.size();
            auto& r = ref[subjId];

            auto edit = std::vector<char>{};
            edit.resize(l, 'M');
            while (e > 0) {
                auto t = rand() % 3;
                if (t == 0) { // substitution
                    auto pos2 = rand() % edit.size();
                    if (edit[pos2] != 'M') continue;
                    edit[pos2] = 'S';
                    --e;
                } else if (t == 1) { // insertion
                    auto pos2 = rand() % edit.size();
                    if (edit[pos2] != 'M') continue;
                    edit[pos2] = 'I';
                    --e;
                } else if (t == 2) { // deletion
                    auto pos2 = rand() % (edit.size()+1);
                    edit.insert(edit.begin() + pos2, 'D');
                    --e;
                }
            }

            auto pos = rand() % (r.size()-l-e);
            auto read = std::vector<uint8_t>{};
            for (auto t : edit) {
                if (t == 'M') {
                    read.push_back(r[pos]);
                    ++pos;
                } else if (t == 'S') {
                    auto c = ((r[pos]-1) + (rand() % 3)) % 4 + 1;
                    read.push_back(c);
                    ++pos;
                } else if (t == 'I') {
                    auto c = (rand() % 4) + 1;
                    read.push_back(c);
                } else if (t == 'D') {
                    ++pos;
                }
            }
            auto gt = std::make_tuple(subjId, pos);
            return std::make_tuple(gt, read);
        };

        auto generateReads = [generateRead](size_t n, std::vector<std::vector<uint8_t>> const& ref, size_t l, size_t e) {
            auto reads = std::vector<std::vector<uint8_t>>{};
            auto expected = std::vector<std::tuple<size_t, size_t, size_t>>{};
            reads.reserve(n);
            for (size_t i{0}; i < n; ++i) {
                auto [gt, read] = generateRead(ref, l, e);
                reads.emplace_back(read);
                expected.emplace_back(i, std::get<0>(gt), std::get<1>(gt));
            }
            return std::make_tuple(expected, reads);
        };

        // generate reference
        auto ref = generateSequences(100, 1'000'000);

        // generate reads
        size_t errors = 2;
        size_t len = 150;
        auto [expected, reads] = generateReads(1000, ref, len, errors);

        auto index = Index{ref, /*samplingRate*/1, /*threadNbr*/1};
        static auto bench = ankerl::nanobench::Bench();
        bench.batch(reads.size())
             .relative(true);

        size_t r_ng12{};
        size_t r_ng21{};
        {
            auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, errors), len);

            bench.run("search ng12", [&]() {
                r_ng12 = 0;
                fmindex_collection::search_ng21::search(index, reads, search_scheme, [&](auto qidx, auto cursor, auto errors) {
                    (void)errors;
                    (void)qidx;
                    r_ng12 += cursor.count();
                    ankerl::nanobench::doNotOptimizeAway(cursor);
                });
            });
        }

        {
            auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, errors), len);

            bench.run("search ng21", [&]() {
                r_ng21 = 0;
                fmindex_collection::search_ng21::search(index, reads, search_scheme, [&](auto qidx, auto cursor, auto errors) {
                    (void)errors;
                    (void)qidx;
                    r_ng21 += cursor.count();
                    ankerl::nanobench::doNotOptimizeAway(cursor);
                });
            });
        }
    }
}
