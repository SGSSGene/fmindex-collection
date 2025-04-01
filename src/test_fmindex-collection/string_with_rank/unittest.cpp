// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include "utils.h"

TEST_CASE("check if rank on the symbol vectors is working, all sizes", "[string_with_rank][all_sizes]") {
    auto testSigma = []<size_t Sigma>() {
        INFO("Sigma " << Sigma);
        call_with_templates<
            ALLRANKVECTORS(Sigma)>([&]<typename Vector>() {
            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto text = generateText<0, Vector::Sigma>(1000000);

            auto vec = Vector{std::span{text}};
            REQUIRE(vec.size() == text.size());
            {
                for (size_t i{0}; i < text.size(); ++i) {
                    INFO(i);
                    CHECK(vec.symbol(i) == text.at(i));
                }
            }
            {
                auto countRank = [&](size_t idx, uint8_t sym) {
                    size_t acc{};
                    for (size_t i{0}; i < idx; ++i) {
                        acc = acc + (text[i] == sym);
                    }
                    return acc;
                };
                for (size_t i{0}; i <= text.size(); ++i) {
                    INFO(i);
                    for (size_t symb{}; symb < Vector::Sigma; ++symb) {
                        CHECK(vec.rank(i, symb) == countRank(i, symb));
                    }
                }
            }
        });
    };

    SECTION("test different sizes of alphabets") {
        testSigma.operator()<4>();
        testSigma.operator()<5>();
        testSigma.operator()<6>();
        testSigma.operator()<16>();
        testSigma.operator()<255>();
    }

}
TEST_CASE("hand counted, test with 255 alphabet", "[string_with_rank][255][small]") {

    auto text = std::vector<uint8_t>{'H', 'a', 'l', 'l', 'o', ' ', 'W', 'e', 'l', 't'};

    SECTION("checks") {
        call_with_templates<
            ALLRANKVECTORS(255)>([&]<typename Vector>() {
            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto vec = Vector{std::span{text}};
            // checking construction size
            REQUIRE(vec.size() == text.size());

            // checking symbols
            for (size_t i{0}; i < text.size(); ++i) {
                INFO(i);
                CHECK(vec.symbol(i) == text.at(i));
            }
            // test complete vector on symbol()"
            CHECK(vec.symbol( 0) == 'H');
            CHECK(vec.symbol( 1) == 'a');
            CHECK(vec.symbol( 2) == 'l');
            CHECK(vec.symbol( 3) == 'l');
            CHECK(vec.symbol( 4) == 'o');
            CHECK(vec.symbol( 5) == ' ');
            CHECK(vec.symbol( 6) == 'W');
            CHECK(vec.symbol( 7) == 'e');
            CHECK(vec.symbol( 8) == 'l');
            CHECK(vec.symbol( 9) == 't');

            // test complete vector on rank()
            CHECK(vec.rank( 0, ' ') == 0);
            CHECK(vec.rank( 1, ' ') == 0);
            CHECK(vec.rank( 2, ' ') == 0);
            CHECK(vec.rank( 3, ' ') == 0);
            CHECK(vec.rank( 4, ' ') == 0);
            CHECK(vec.rank( 5, ' ') == 0);
            CHECK(vec.rank( 6, ' ') == 1);
            CHECK(vec.rank( 7, ' ') == 1);
            CHECK(vec.rank( 8, ' ') == 1);
            CHECK(vec.rank( 9, ' ') == 1);
            CHECK(vec.rank(10, ' ') == 1);

            CHECK(vec.rank( 0, 'H') == 0);
            CHECK(vec.rank( 1, 'H') == 1);
            CHECK(vec.rank( 2, 'H') == 1);
            CHECK(vec.rank( 3, 'H') == 1);
            CHECK(vec.rank( 4, 'H') == 1);
            CHECK(vec.rank( 5, 'H') == 1);
            CHECK(vec.rank( 6, 'H') == 1);
            CHECK(vec.rank( 7, 'H') == 1);
            CHECK(vec.rank( 8, 'H') == 1);
            CHECK(vec.rank( 9, 'H') == 1);
            CHECK(vec.rank(10, 'H') == 1);

            CHECK(vec.rank( 0, 'W') == 0);
            CHECK(vec.rank( 1, 'W') == 0);
            CHECK(vec.rank( 2, 'W') == 0);
            CHECK(vec.rank( 3, 'W') == 0);
            CHECK(vec.rank( 4, 'W') == 0);
            CHECK(vec.rank( 5, 'W') == 0);
            CHECK(vec.rank( 6, 'W') == 0);
            CHECK(vec.rank( 7, 'W') == 1);
            CHECK(vec.rank( 8, 'W') == 1);
            CHECK(vec.rank( 9, 'W') == 1);
            CHECK(vec.rank(10, 'W') == 1);

            CHECK(vec.rank( 0, 'a') == 0);
            CHECK(vec.rank( 1, 'a') == 0);
            CHECK(vec.rank( 2, 'a') == 1);
            CHECK(vec.rank( 3, 'a') == 1);
            CHECK(vec.rank( 4, 'a') == 1);
            CHECK(vec.rank( 5, 'a') == 1);
            CHECK(vec.rank( 6, 'a') == 1);
            CHECK(vec.rank( 7, 'a') == 1);
            CHECK(vec.rank( 8, 'a') == 1);
            CHECK(vec.rank( 9, 'a') == 1);
            CHECK(vec.rank(10, 'a') == 1);

            CHECK(vec.rank( 0, 'e') == 0);
            CHECK(vec.rank( 1, 'e') == 0);
            CHECK(vec.rank( 2, 'e') == 0);
            CHECK(vec.rank( 3, 'e') == 0);
            CHECK(vec.rank( 4, 'e') == 0);
            CHECK(vec.rank( 5, 'e') == 0);
            CHECK(vec.rank( 6, 'e') == 0);
            CHECK(vec.rank( 7, 'e') == 0);
            CHECK(vec.rank( 8, 'e') == 1);
            CHECK(vec.rank( 9, 'e') == 1);
            CHECK(vec.rank(10, 'e') == 1);

            CHECK(vec.rank( 0, 'l') == 0);
            CHECK(vec.rank( 1, 'l') == 0);
            CHECK(vec.rank( 2, 'l') == 0);
            CHECK(vec.rank( 3, 'l') == 1);
            CHECK(vec.rank( 4, 'l') == 2);
            CHECK(vec.rank( 5, 'l') == 2);
            CHECK(vec.rank( 6, 'l') == 2);
            CHECK(vec.rank( 7, 'l') == 2);
            CHECK(vec.rank( 8, 'l') == 2);
            CHECK(vec.rank( 9, 'l') == 3);
            CHECK(vec.rank(10, 'l') == 3);

            CHECK(vec.rank( 0, 'o') == 0);
            CHECK(vec.rank( 1, 'o') == 0);
            CHECK(vec.rank( 2, 'o') == 0);
            CHECK(vec.rank( 3, 'o') == 0);
            CHECK(vec.rank( 4, 'o') == 0);
            CHECK(vec.rank( 5, 'o') == 1);
            CHECK(vec.rank( 6, 'o') == 1);
            CHECK(vec.rank( 7, 'o') == 1);
            CHECK(vec.rank( 8, 'o') == 1);
            CHECK(vec.rank( 9, 'o') == 1);
            CHECK(vec.rank(10, 'o') == 1);

            CHECK(vec.rank( 0, 't') == 0);
            CHECK(vec.rank( 1, 't') == 0);
            CHECK(vec.rank( 2, 't') == 0);
            CHECK(vec.rank( 3, 't') == 0);
            CHECK(vec.rank( 4, 't') == 0);
            CHECK(vec.rank( 5, 't') == 0);
            CHECK(vec.rank( 6, 't') == 0);
            CHECK(vec.rank( 7, 't') == 0);
            CHECK(vec.rank( 8, 't') == 0);
            CHECK(vec.rank( 9, 't') == 0);
            CHECK(vec.rank(10, 't') == 1);

            // check all other characters have rank 0
            auto ignore = std::unordered_set<size_t>{' ', 'H', 'W', 'a', 'e', 'l', 'o', 't'};
            for (size_t s{0}; s < Vector::Sigma; ++s) {
                if (ignore.contains(s)) continue;
                for (size_t i{0}; i < 11; ++i) {
                    CHECK(vec.rank(i, s) == 0);
                }
            }

            // test complete vec 'H' for prefix_rank()
            CHECK(vec.prefix_rank( 0, ' ') == 0);
            CHECK(vec.prefix_rank( 1, ' ') == 0);
            CHECK(vec.prefix_rank( 2, ' ') == 0);
            CHECK(vec.prefix_rank( 3, ' ') == 0);
            CHECK(vec.prefix_rank( 4, ' ') == 0);
            CHECK(vec.prefix_rank( 5, ' ') == 0);
            CHECK(vec.prefix_rank( 6, ' ') == 1);
            CHECK(vec.prefix_rank( 7, ' ') == 1);
            CHECK(vec.prefix_rank( 8, ' ') == 1);
            CHECK(vec.prefix_rank( 9, ' ') == 1);
            CHECK(vec.prefix_rank(10, ' ') == 1);

            CHECK(vec.prefix_rank( 0, 'H') == 0);
            CHECK(vec.prefix_rank( 1, 'H') == 1);
            CHECK(vec.prefix_rank( 2, 'H') == 1);
            CHECK(vec.prefix_rank( 3, 'H') == 1);
            CHECK(vec.prefix_rank( 4, 'H') == 1);
            CHECK(vec.prefix_rank( 5, 'H') == 1);
            CHECK(vec.prefix_rank( 6, 'H') == 2);
            CHECK(vec.prefix_rank( 7, 'H') == 2);
            CHECK(vec.prefix_rank( 8, 'H') == 2);
            CHECK(vec.prefix_rank( 9, 'H') == 2);
            CHECK(vec.prefix_rank(10, 'H') == 2);

            CHECK(vec.prefix_rank( 0, 'W') == 0);
            CHECK(vec.prefix_rank( 1, 'W') == 1);
            CHECK(vec.prefix_rank( 2, 'W') == 1);
            CHECK(vec.prefix_rank( 3, 'W') == 1);
            CHECK(vec.prefix_rank( 4, 'W') == 1);
            CHECK(vec.prefix_rank( 5, 'W') == 1);
            CHECK(vec.prefix_rank( 6, 'W') == 2);
            CHECK(vec.prefix_rank( 7, 'W') == 3);
            CHECK(vec.prefix_rank( 8, 'W') == 3);
            CHECK(vec.prefix_rank( 9, 'W') == 3);
            CHECK(vec.prefix_rank(10, 'W') == 3);

            CHECK(vec.prefix_rank( 0, 'a') == 0);
            CHECK(vec.prefix_rank( 1, 'a') == 1);
            CHECK(vec.prefix_rank( 2, 'a') == 2);
            CHECK(vec.prefix_rank( 3, 'a') == 2);
            CHECK(vec.prefix_rank( 4, 'a') == 2);
            CHECK(vec.prefix_rank( 5, 'a') == 2);
            CHECK(vec.prefix_rank( 6, 'a') == 3);
            CHECK(vec.prefix_rank( 7, 'a') == 4);
            CHECK(vec.prefix_rank( 8, 'a') == 4);
            CHECK(vec.prefix_rank( 9, 'a') == 4);
            CHECK(vec.prefix_rank(10, 'a') == 4);

            CHECK(vec.prefix_rank( 0, 'e') == 0);
            CHECK(vec.prefix_rank( 1, 'e') == 1);
            CHECK(vec.prefix_rank( 2, 'e') == 2);
            CHECK(vec.prefix_rank( 3, 'e') == 2);
            CHECK(vec.prefix_rank( 4, 'e') == 2);
            CHECK(vec.prefix_rank( 5, 'e') == 2);
            CHECK(vec.prefix_rank( 6, 'e') == 3);
            CHECK(vec.prefix_rank( 7, 'e') == 4);
            CHECK(vec.prefix_rank( 8, 'e') == 5);
            CHECK(vec.prefix_rank( 9, 'e') == 5);
            CHECK(vec.prefix_rank(10, 'e') == 5);

            CHECK(vec.prefix_rank( 0, 'l') == 0);
            CHECK(vec.prefix_rank( 1, 'l') == 1);
            CHECK(vec.prefix_rank( 2, 'l') == 2);
            CHECK(vec.prefix_rank( 3, 'l') == 3);
            CHECK(vec.prefix_rank( 4, 'l') == 4);
            CHECK(vec.prefix_rank( 5, 'l') == 4);
            CHECK(vec.prefix_rank( 6, 'l') == 5);
            CHECK(vec.prefix_rank( 7, 'l') == 6);
            CHECK(vec.prefix_rank( 8, 'l') == 7);
            CHECK(vec.prefix_rank( 9, 'l') == 8);
            CHECK(vec.prefix_rank(10, 'l') == 8);

            CHECK(vec.prefix_rank( 0, 'o') == 0);
            CHECK(vec.prefix_rank( 1, 'o') == 1);
            CHECK(vec.prefix_rank( 2, 'o') == 2);
            CHECK(vec.prefix_rank( 3, 'o') == 3);
            CHECK(vec.prefix_rank( 4, 'o') == 4);
            CHECK(vec.prefix_rank( 5, 'o') == 5);
            CHECK(vec.prefix_rank( 6, 'o') == 6);
            CHECK(vec.prefix_rank( 7, 'o') == 7);
            CHECK(vec.prefix_rank( 8, 'o') == 8);
            CHECK(vec.prefix_rank( 9, 'o') == 9);
            CHECK(vec.prefix_rank(10, 'o') == 9);

            CHECK(vec.prefix_rank( 0, 't') ==  0);
            CHECK(vec.prefix_rank( 1, 't') ==  1);
            CHECK(vec.prefix_rank( 2, 't') ==  2);
            CHECK(vec.prefix_rank( 3, 't') ==  3);
            CHECK(vec.prefix_rank( 4, 't') ==  4);
            CHECK(vec.prefix_rank( 5, 't') ==  5);
            CHECK(vec.prefix_rank( 6, 't') ==  6);
            CHECK(vec.prefix_rank( 7, 't') ==  7);
            CHECK(vec.prefix_rank( 8, 't') ==  8);
            CHECK(vec.prefix_rank( 9, 't') ==  9);
            CHECK(vec.prefix_rank(10, 't') == 10);

            // check all_ranks() is equal to prefix_rank() and rank()
            for (size_t idx{0}; idx < vec.size(); ++idx) {
                auto [rank, prefix] = vec.all_ranks_and_prefix_ranks(idx);
                auto rank2 = vec.all_ranks(idx);
                for (size_t symb{1}; symb < Vector::Sigma; ++symb) {
                    INFO(idx);
                    INFO(symb);
                    CHECK(rank[symb] == vec.rank(idx, symb));
                    CHECK(rank2[symb] == vec.rank(idx, symb));
                    CHECK(prefix[symb] == vec.prefix_rank(idx, symb));
                }
            }
        });
    }
}

TEST_CASE("check symbol vectors construction on text longer than 255 characters", "[string_with_rank][255][large]") {
    auto text = std::vector<uint8_t>{'H', 'a', 'l', 'l', 'o', ' ', 'W', 'e', 'l', 't',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     254, 254, 254, 254, 254, 254, 254, 254, 254, 254,
                                    };

    SECTION("checks") {
        call_with_templates<
            ALLRANKVECTORS(255)>([&]<typename Vector>() {

            auto vector = Vector{text};

            REQUIRE(vector.size() == text.size());

            // check that symbol() call works
            for (size_t i{0}; i < text.size(); ++i) {
                INFO(i);
                CHECK(vector.symbol(i) == text.at(i));
            }

            auto countRank = [&](size_t idx, uint8_t sym) {
                size_t acc{};
                for (size_t i{0}; i < idx; ++i) {
                    acc = acc + (text[i] == sym);
                }
                return acc;
            };
            auto countPrefixRank = [&](size_t idx, uint8_t sym) {
                size_t acc{};
                for (size_t i{0}; i < idx; ++i) {
                    acc = acc + (text[i] <= sym);
                }
                return acc;
            };


            // check all_ranks() is equal to prefix_rank() and rank()
            for (size_t idx{0}; idx <= vector.size(); ++idx) {
                auto [rank, prefix] = vector.all_ranks_and_prefix_ranks(idx);
                auto rank2 = vector.all_ranks(idx);
                for (size_t symb{1}; symb < Vector::Sigma; ++symb) {
                    INFO(idx);
                    INFO(symb);
                    CHECK(countRank(idx, symb) == vector.rank(idx, symb));
                    CHECK(countPrefixRank(idx, symb) == vector.prefix_rank(idx, symb));
                    vector.rank(idx, symb);
                    vector.prefix_rank(idx, symb);
                    CHECK(rank[symb] == vector.rank(idx, symb));
                    CHECK(rank2[symb] == vector.rank(idx, symb));
                    CHECK(countRank(idx, symb) == rank[symb]);
                    CHECK(countPrefixRank(idx, symb) == prefix[symb]);
                    CHECK(prefix[symb] == vector.prefix_rank(idx, symb));
                }
            }
        });
    }
}

#if 0
// specific benchmarks, that load from a file
static auto benchs_6 = Benchs{};
static auto benchs_6_text = std::vector<uint8_t>{};
static auto benchSize_bwt = BenchSize{};
TEMPLATE_TEST_CASE("benchmark vectors c'tor operation, on human dna5 data", "[RankVector][bwt][!benchmark][6][time][ctor][.]", ALLRANKVECTORS(6)) {
    using Vector = TestType;
    if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<6>>) {
        return;
    }
    auto& [bench_rank, bench_prefix_rank, bench_all_ranks, bench_all_prefix_ranks, bench_symbol, bench_ctor] = benchs_6;


    auto vector_name = getName<Vector>();
    INFO(vector_name);

    SECTION("benchmarking") {
        auto rng = ankerl::nanobench::Rng{};

        // generates string with values between 1-4
        auto& text = benchs_6_text;
        if (text.empty()) {
            auto fs = std::ifstream{"bwt.bwt"};
            std::string line;
            std::getline(fs, line);
            fs.close();
            text.resize(line.size());
            std::memcpy(text.data(), line.data(), line.size());
            for (auto& c : text) {
                switch(c) {
                case '$': c = 0; break;
                case 'A': c = 1; break;
                case 'C': c = 2; break;
                case 'G': c = 3; break;
                case 'T': c = 4; break;
                case 'N': c = 5; break;
                default:
                    throw std::runtime_error{"unknown char " + std::to_string((int)c)};

                }
            }
        }

        auto vec = Vector{text};


        size_t minEpochIterations = 2'000'000;
        minEpochIterations = 1;
        bench_ctor.minEpochIterations(minEpochIterations).batch(text.size()).run(vector_name, [&]() {
            auto vec = Vector{text};
            ankerl::nanobench::doNotOptimizeAway(const_cast<Vector const&>(vec));
        });
    }
}

TEMPLATE_TEST_CASE("benchmark vectors symbol() and rank() operations, on human dna5 data", "[RankVector][bwt][!benchmark][6][time][.]", ALLRANKVECTORS(6)) {
    using Vector = TestType;
    if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<6>>) {
        return;
    }
    auto& [bench_rank, bench_prefix_rank, bench_all_ranks, bench_all_prefix_ranks, bench_symbol, bench_ctor] = benchs_6;


    auto vector_name = getName<Vector>();
    INFO(vector_name);

    SECTION("benchmarking") {
        auto rng = ankerl::nanobench::Rng{};

        // generates string with values between 1-4
        auto& text = benchs_6_text;
        if (text.empty()) {
            auto fs = std::ifstream{"bwt.bwt"};
            std::string line;
            std::getline(fs, line);
            fs.close();
            text.resize(line.size());
            std::memcpy(text.data(), line.data(), line.size());
            for (auto& c : text) {
                switch(c) {
                case '$': c = 0; break;
                case 'A': c = 1; break;
                case 'C': c = 2; break;
                case 'G': c = 3; break;
                case 'T': c = 4; break;
                case 'N': c = 5; break;
                default:
                    throw std::runtime_error{"unknown char " + std::to_string((int)c)};

                }
            }
        }

        auto vec = Vector{text};
        REQUIRE(text.size() > 0);


        size_t minEpochIterations = 2'000'000;
        minEpochIterations = 1;
        bench_symbol.minEpochIterations(minEpochIterations).run(vector_name, [&]() {
            auto v = vec.symbol(rng.bounded(text.size()));
            ankerl::nanobench::doNotOptimizeAway(v);
        });

        bench_rank.minEpochIterations(minEpochIterations).run(vector_name, [&]() {
            auto v = vec.rank(rng.bounded(text.size()), rng.bounded(5)+1);
            ankerl::nanobench::doNotOptimizeAway(v);
        });

        bench_prefix_rank.minEpochIterations(minEpochIterations).run(vector_name, [&]() {
            auto v = vec.prefix_rank(rng.bounded(text.size()), rng.bounded(5)+1);
            ankerl::nanobench::doNotOptimizeAway(v);
        });

        bench_all_ranks.minEpochIterations(minEpochIterations).run(vector_name, [&]() {
            auto v = vec.all_ranks(rng.bounded(text.size()));
            ankerl::nanobench::doNotOptimizeAway(v);
        });

        bench_all_prefix_ranks.minEpochIterations(minEpochIterations).run(vector_name, [&]() {
            auto v = vec.all_ranks_and_prefix_ranks(rng.bounded(text.size()));
            ankerl::nanobench::doNotOptimizeAway(v);
        });
    }
}

TEMPLATE_TEST_CASE("benchmark vectors size, on human dna5 data", "[RankVector][bwt][!benchmark][6][size][.]", ALLRANKVECTORS(6)) {
    using Vector = TestType;
    if constexpr (std::same_as<Vector, fmindex_collection::rankvector::Naive<6>>) {
        return;
    }

    auto vector_name = getName<Vector>();
    INFO(vector_name);

    SECTION("benchmarking") {
        auto rng = ankerl::nanobench::Rng{};

        // generates string with values between 1-4
        auto& text = benchs_6_text;
        if (text.empty()) {
            auto fs = std::ifstream{"bwt.bwt"};
            std::string line;
            std::getline(fs, line);
            fs.close();
            text.resize(line.size());
            std::memcpy(text.data(), line.data(), line.size());
            for (auto& c : text) {
                switch(c) {
                case '$': c = 0; break;
                case 'A': c = 1; break;
                case 'C': c = 2; break;
                case 'G': c = 3; break;
                case 'T': c = 4; break;
                case 'N': c = 5; break;
                default:
                    throw std::runtime_error{"unknown char " + std::to_string((int)c)};

                }
            }
        }

        auto vec = Vector{text};

        auto ofs     = std::stringstream{};
        auto archive = cereal::BinaryOutputArchive{ofs};
        archive(vec);
        auto s = ofs.str().size();
        benchSize_bwt.addEntry({
            .name = vector_name,
            .size = s,
            .text_size = text.size(),
            .bits_per_char = (s*8)/double(text.size())
        });
    }
}
#endif
