// SPDX-FileCopyrightText: 2023 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2023 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <catch2/catch_all.hpp>
#include <fmindex-collection/fmindex/BiFMIndex.h>
#include <fmindex-collection/locate.h>
#include <fmindex-collection/search/all.h>
#include <fmindex-collection/search_scheme/generator/all.h>
#include <fmindex-collection/search_scheme/expand.h>
#include <fmindex-collection/string/all.h>
#include <nanobench.h>

TEST_CASE("check searches with errors", "[searches][errors]") {
    using Index = fmc::BiFMIndex<256>;

    auto input  = std::vector<std::vector<uint8_t>>{{'A', 'A', 'A', 'C', 'A', 'A', 'A', 'B', 'A', 'A', 'A'},
                                                    {'A', 'A', 'A', 'B', 'A', 'A', 'A', 'C', 'A', 'A', 'A'}};

    auto index = Index{input, /*samplingRate*/1, /*threadNbr*/1};

    auto queries = std::vector<std::vector<uint8_t>> {std::vector<uint8_t>{'C', 'C'}, std::vector<uint8_t>{'B', 'B'}};
    SECTION("backtracking, single search") {
        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        for (size_t qidx{0}; qidx != queries.size(); ++qidx) {
            auto const& query = queries[qidx];
            fmc::search_backtracking::search(index, query, 1, [&](auto cursor, auto errors) {
                (void)errors;
                for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                    auto [sid, spos] = entry;
                    results.emplace_back(qidx, sid, spos+offset);
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
        fmc::search_backtracking::search(index, queries, 1, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto buffer1 = std::vector<std::pair<fmc::select_cursor_t<Index>, size_t>>{};
        auto buffer2 = std::vector<std::pair<fmc::select_cursor_t<Index>, size_t>>{};

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        for (size_t qidx{0}; qidx < queries.size(); ++qidx) {
            auto const& query = queries[qidx];
            fmc::search_backtracking_with_buffers::search(index, query, 1, buffer1, buffer2, [&](auto cursor, auto errors) {
                (void)errors;
                for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                    auto [sid, spos] = entry;

                    results.emplace_back(qidx, sid, spos+offset);
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
        fmc::search_no_errors::search(index, queries, [&](auto qidx, auto cursor) {
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
        };
        CHECK(results == expected);
    }

    SECTION("search one error, all search") {
        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_one_error::search(index, queries, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_pseudo::search</*EditDistance=*/true>(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        for (size_t qidx{0}; qidx < queries.size(); ++qidx) {
            auto const& query = queries[qidx];
            fmc::search_pseudo::search</*EditDistance=*/true>(index, query, search_scheme, [&](auto cursor, auto errors) {
                (void)errors;
                for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                    auto [sid, spos] = entry;
                    results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_pseudo::search</*EditDistance=*/false>(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng12::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng14::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng15::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::backtracking(2, 0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng16::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng17::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::backtracking(2, 0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng20::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng21::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng21::search_n(index, queries, search_scheme, 3, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme_e0 = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 0), queries[0].size());
        auto search_scheme_e1 = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(1, 1), queries[0].size());
        auto search_scheme_e2 = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(2, 2), queries[0].size());
        auto search_schemes = std::vector{search_scheme_e0, search_scheme_e1, search_scheme_e2};

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng21::search_best(index, queries, search_schemes, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme_e0 = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 0), queries[0].size());
        auto search_scheme_e1 = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(1, 1), queries[0].size());
        auto search_schemes = std::vector{search_scheme_e0, search_scheme_e1};

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng21::search_best_n(index, queries, search_schemes, 3, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng21V2::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng21V3::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng21V4::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng21V5::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng21V6::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng21V6::search_n(index, queries, search_scheme, 3, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme_e0 = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 0), queries[0].size());
        auto search_scheme_e1 = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(1, 1), queries[0].size());
        auto search_scheme_e2 = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(2, 2), queries[0].size());
        auto search_schemes = std::vector{search_scheme_e0, search_scheme_e1, search_scheme_e2};

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng21V6::search_best(index, queries, search_schemes, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme_e0 = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 0), queries[0].size());
        auto search_scheme_e1 = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(1, 1), queries[0].size());
        auto search_schemes = std::vector{search_scheme_e0, search_scheme_e1};

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng21V6::search_best_n(index, queries, search_schemes, 3, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng21V7::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng21V7::search_n(index, queries, search_scheme, 3, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 2), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng21V7::search_best(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 2), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng21V7::search_best_n(index, queries, search_scheme, 3, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng22::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors, auto const& action) {
            (void)errors;
            (void)action;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
    SECTION("search ng24, all search") {
        auto input  = std::vector<std::vector<uint8_t>>{{'A', 'A', 'A', 'C', 'A', 'A', 'A', 'B', 'A', 'A', 'A'},
                                                        {'A', 'A', 'A', 'B', 'A', 'A', 'A', 'C', 'A', 'A', 'A'}};

        auto index = Index{input, /*samplingRate*/1, /*threadNbr*/1};

        auto queries = std::vector<std::vector<uint8_t>> {std::vector<uint8_t>{'C', 'D'}, std::vector<uint8_t>{'D', 'B'}};

        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng24::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 1, 7},
            {1, 0, 7},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }
    SECTION("search ng24, all search_n") {
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng24::search(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
            }
        }, 3);

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

    SECTION("search ng24, hamming distance, all search") {
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng24::search</*EditDistance=*/false>(index, queries, search_scheme, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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

    SECTION("search ng24, all search, no search scheme") {
        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng24::search(index, queries, /*.maxErrors=*/1, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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

    SECTION("search ng24, all search_n, no search scheme") {
        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng24::search(index, queries, /*.maxErrors=*/1, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
            }
        }, /*.n=*/3);

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

    SECTION("search ng24, hamming distance, all search, no search scheme") {
        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng24::search</*EditDistance=*/false>(index, queries, /*.maxErrors=*/1, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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

    SECTION("search ng25, all search") {
        auto input  = std::vector<std::vector<uint8_t>>{{'A', 'A', 'A', 'C', 'A', 'A', 'A', 'B', 'A', 'A', 'A'},
                                                        {'A', 'A', 'A', 'B', 'A', 'A', 'A', 'C', 'A', 'A', 'A'}};

        auto index = Index{input, /*samplingRate*/1, /*threadNbr*/1};

        auto queries = std::vector<std::vector<uint8_t>> {std::vector<uint8_t>{'C', 'D'}, std::vector<uint8_t>{'D', 'B'}};

        auto search_scheme = fmc::search_scheme::generator::pigeon_opt(0, 1);
        auto partition     = fmc::search_scheme::createUniformPartition(search_scheme, queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng25::search(index, queries, search_scheme, partition, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 1, 7},
            {1, 0, 7},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }
    SECTION("search ng25, all search_n") {
        auto search_scheme = fmc::search_scheme::generator::pigeon_opt(0, 1);
        auto partition     = fmc::search_scheme::createUniformPartition(search_scheme, queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng25::search(index, queries, search_scheme, partition, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
            }
        }, 3);

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

    SECTION("search ng25, hamming distance, all search") {
        auto search_scheme = fmc::search_scheme::generator::pigeon_opt(0, 1);
        auto partition     = fmc::search_scheme::createUniformPartition(search_scheme, queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng25::search</*EditDistance=*/false>(index, queries, search_scheme, partition, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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

    SECTION("search ng26, all search") {
        auto input  = std::vector<std::vector<uint8_t>>{{'A', 'A', 'A', 'C', 'A', 'A', 'A', 'B', 'A', 'A', 'A'},
                                                        {'A', 'A', 'A', 'B', 'A', 'A', 'A', 'C', 'A', 'A', 'A'}};

        auto index = Index{input, /*samplingRate*/1, /*threadNbr*/1};

        auto queries = std::vector<std::vector<uint8_t>> {std::vector<uint8_t>{'C', 'D'}, std::vector<uint8_t>{'D', 'B'}};

        auto search_scheme = fmc::search_scheme::generator::pigeon_opt(0, 1);
        auto partition     = fmc::search_scheme::createUniformPartition(search_scheme, queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng26::search(index, queries, search_scheme, partition, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 1, 7},
            {1, 0, 7},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }
    #if 0
    SECTION("search, all search, no search scheme") {
        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search</*EditDistance=*/true>(index, queries, /*.maxErrors=*/1, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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
    #endif

    SECTION("search ng26, all search_n") {
        auto search_scheme = fmc::search_scheme::generator::pigeon_opt(0, 1);
        auto partition     = fmc::search_scheme::createUniformPartition(search_scheme, queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng26::search(index, queries, search_scheme, partition, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
            }
        }, 3);

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

    SECTION("search ng26, hamming distance, all search") {
        auto search_scheme = fmc::search_scheme::generator::pigeon_opt(0, 1);
        auto partition     = fmc::search_scheme::createUniformPartition(search_scheme, queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng25::search</*EditDistance=*/false>(index, queries, search_scheme, partition, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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



    SECTION("search double index, all search") {
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());

        auto queryIndex = Index{queries, /*samplingRate*/1, /*threadNbr*/1};

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_double_index::search<true>(index, queryIndex, search_scheme, /*.threshold=*/1, /*.optMode=*/7, [&](auto cursor, auto qcursor, auto errors) {
            (void)errors;
            auto qidxs = std::vector<size_t>{};
            for (auto [entry, offset] : fmc::LocateLinear{queryIndex, qcursor}) {
                auto [sid, spos] = entry;
                qidxs.push_back(sid);
            }

            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                for (auto qidx : qidxs) {
                    results.emplace_back(qidx, sid, spos+offset);
                }
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 1, 7},
            {1, 0, 7},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("search double index 2, all search") {
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};

        using QIndex = fmc::LinearFMIndex<Index::Sigma>;
        auto ordered_queries = queries;

        //!TODO this function is not really ready (also doesn't deliver on the speed promises)
        for (auto const& search : search_scheme) {
            // bring queries parts into correct order
            for (size_t qid{0}; qid < queries.size(); ++qid) {
                for (size_t i{0}; i < search.pi.size(); ++i) {
                    ordered_queries[qid][search.pi.size() - i - 1] = queries[qid][search.pi[i]];
                }
            }
            auto queryIndex = QIndex{ordered_queries};

            fmc::search_double_index2::search<true>(index, queryIndex, search, /*.threshold=*/1, /*.optMode=*/7, [&](auto cursor, auto qcursor, auto errors) {
                (void)errors;
                auto qidxs = std::vector<size_t>{};
                for (auto row : qcursor) {
                    qidxs.push_back(queryIndex.locate(row));
                }

                for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                    auto [sid, spos] = entry;

                    for (auto qidx : qidxs) {
                        results.emplace_back(qidx, sid, spos+offset);
                    }
                }
            });
        }

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 1, 7},
            {1, 0, 7},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("search ng24sm, all search") {
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());
        auto scoringMatrix = fmc::search_ng24sm::ScoringMatrix<Index::Sigma>{/*ambiguous=*/0};

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng24sm::search(index, queries, search_scheme, scoringMatrix, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 1, 7},
            {1, 0, 7},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("search ng24sm, hamming distance, all search") {
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());
        auto scoringMatrix = fmc::search_ng24sm::ScoringMatrix<Index::Sigma>{/*ambiguous=*/0};

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng24sm::search</*EditDistance=*/false>(index, queries, search_scheme, scoringMatrix, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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


    SECTION("search ng24sm, all search_n") {
        auto search_scheme = fmc::search_scheme::expand(fmc::search_scheme::generator::pigeon_opt(0, 1), queries[0].size());
        auto scoringMatrix = fmc::search_ng24sm::ScoringMatrix<Index::Sigma>{/*ambiguous=*/0};

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng24sm::search_n(index, queries, search_scheme, /*n=*/3, scoringMatrix, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 1, 7},
            {1, 0, 7},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("search ng24sm, all search, no search scheme") {
        auto scoringMatrix = fmc::search_ng24sm::ScoringMatrix<Index::Sigma>{/*ambiguous=*/0};

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng24sm::search(index, queries, /*.maxErrors=*/1, scoringMatrix, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 1, 7},
            {1, 0, 7},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("search ng24sm, hamming distance, all search, no search scheme") {
        auto scoringMatrix = fmc::search_ng24sm::ScoringMatrix<Index::Sigma>{/*ambiguous=*/0};

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng24sm::search</*EditDistance=*/false>(index, queries, /*.maxErrors=*/1, scoringMatrix, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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


    SECTION("search ng24sm, all search_n, no search scheme") {
        auto scoringMatrix = fmc::search_ng24sm::ScoringMatrix<Index::Sigma>{/*ambiguous=*/0};

        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_ng24sm::search_n(index, queries, /*.maxErrors=*/1, /*.n=*/3, scoringMatrix, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
            }
        });

        std::ranges::sort(results);

        auto expected = std::vector<std::tuple<size_t, size_t, size_t>> {
            {0, 0, 3},
            {0, 1, 7},
            {1, 0, 7},
            {1, 1, 3},
        };
        CHECK(results == expected);
    }

    SECTION("search, all search, no search scheme") {
        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search</*EditDistance=*/true>(index, queries, /*.maxErrors=*/1, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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

    SECTION("search, all search_n, no search scheme") {
        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search_n</*EditDistance=*/true>(index, queries, /*.maxErrors=*/1, /*.n=*/3, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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

    SECTION("search, hamming distance, all search, no search scheme") {
        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::search</*EditDistance=*/false>(index, queries, /*.maxErrors=*/1, [&](auto qidx, auto cursor, auto errors) {
            (void)errors;
            for (auto [entry, offset] : fmc::LocateLinear{index, cursor}) {
                auto [sid, spos] = entry;

                results.emplace_back(qidx, sid, spos+offset);
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

    SECTION("simple search, all search") {
        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::Search {
            .index      = index,
            .queries    = queries,
            .errors     = 1,
            .reportFunc = [&](auto qidx, auto sid, auto spos, auto errors) {
                (void)errors;
                results.emplace_back(qidx, sid, spos);
            },
        }();

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

    SECTION("simple search, all search_n") {
        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::Search {
            .index      = index,
            .queries    = queries,
            .errors     = 1,
            .maxResults = 3,
            .reportFunc = [&](auto qidx, auto sid, auto spos, auto errors) {
                (void)errors;
                results.emplace_back(qidx, sid, spos);
            },
        }();

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

    SECTION("simple search, hamming distance, all search") {
        auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
        fmc::Search {
            .index        = index,
            .queries      = queries,
            .editDistance = false,
            .errors       = 1,
            .reportFunc   = [&](auto qidx, auto sid, auto spos, auto errors) {
                (void)errors;
                results.emplace_back(qidx, sid, spos);
            },
        }();

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
