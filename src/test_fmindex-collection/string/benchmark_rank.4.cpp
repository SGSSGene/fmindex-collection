// SPDX-FileCopyrightText: 2026 Simon Gene Gottlieb
// SPDX-License-Identifier: CC0-1.0
#include "allStrings.h"
#include "utils.h"

#include <fstream>

template <size_t Sigma>
using String = fmc::string::PairedFlattenedBitvectors_512_64k<Sigma>;

template <size_t Sigma>
using StringPartialSymb = fmc::string::PairedFlattenedBitvectorsPartialSymb_512_64k<Sigma>;

template <size_t Sigma>
using String2 = fmc::string::InterleavedBitvectorPrefix16<Sigma>;


TEST_CASE("benchmark string rank/prefix_rank with _limit", "[string][!benchmark][rank][prefix_rank]") {

    static auto const& text16 = generateText<0, 16>();
    static auto const& text4  = [&]() {
        static std::string ret;
        ret.reserve(text16.size());
        for (auto v : text16) {
            ret.push_back(v >> 2);
        }
        return ret;
    }();



    auto s16 = String<16>(text16);
    auto s4  = String<4>(text4);
    auto s16_b = StringPartialSymb<16>(text16);
    auto s4_b  = StringPartialSymb<4>(text4);
    auto s16_i = String2<16>(text16);
    auto s4_i  = String2<16>(text16);

#if 1
    SECTION("benchmarking rank") {
        auto rng = ankerl::nanobench::Rng{};

        auto bench = ankerl::nanobench::Bench{};
        bench.title("rank")
             .relative(true)
             .batch(100);

        bench.run("rank-4", [&]() {
            auto symb = rng.bounded(4);
            size_t a{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size());
                a += s4.rank(pos, symb);
            }
            ankerl::nanobench::doNotOptimizeAway(a);

        });
        bench.run("rank_limit<2>-4", [&]() {
            auto symb = rng.bounded(4);
            size_t a{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size());
                a += s4.rank_limit<2>(pos, symb);
            }
            ankerl::nanobench::doNotOptimizeAway(a);
        });
        bench.run("rank-16", [&]() {
            auto symb = rng.bounded(16);
            size_t a{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text16.size());
                a += s16.rank(pos, symb);
            }
            ankerl::nanobench::doNotOptimizeAway(a);
        });
        bench.run("rank_limit<2>-16", [&]() {
            auto symb = rng.bounded(4);
            size_t a{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text16.size());
                a += s16.rank_limit<2>(pos, symb);
            }
            ankerl::nanobench::doNotOptimizeAway(a);
        });

        bench.run("b rank-4", [&]() {
            auto symb = rng.bounded(4);
            size_t a{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size());
                a += s4_b.rank(pos, symb);
            }
            ankerl::nanobench::doNotOptimizeAway(a);

        });
        bench.run("b rank_limit<2>-4", [&]() {
            auto symb = rng.bounded(4);
            size_t a{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size());
                a += s4_b.rank_limit(pos, symb);
            }
            ankerl::nanobench::doNotOptimizeAway(a);
        });
        bench.run("b rank-16", [&]() {
            auto symb = rng.bounded(16);
            size_t a{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text16.size());
                a += s16_b.rank(pos, symb);
            }
            ankerl::nanobench::doNotOptimizeAway(a);
        });
        bench.run("b rank_limit<2>-16", [&]() {
            auto symb = rng.bounded(4);
            size_t a{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text16.size());
                a += s16_b.rank_limit(pos, symb);
            }
            ankerl::nanobench::doNotOptimizeAway(a);
        });

        bench.run("interleaved rank-4", [&]() {
            auto symb = rng.bounded(4);
            size_t a{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size());
                a += s4_i.rank(pos, symb);
            }
            ankerl::nanobench::doNotOptimizeAway(a);

        });
        bench.run("interleaved rank_limit<2>-4", [&]() {
            auto symb = rng.bounded(4);
            size_t a{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size());
                a += s4_i.rank_limit<2>(pos, symb);
            }
            ankerl::nanobench::doNotOptimizeAway(a);

        });

        bench.run("interleaved rank-16", [&]() {
            auto symb = rng.bounded(16);
            size_t a{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text16.size());
                a += s16_i.rank(pos, symb);
            }
            ankerl::nanobench::doNotOptimizeAway(a);
        });
        bench.run("interleaved rank_limit<2>-16", [&]() {
            auto symb = rng.bounded(4);
            size_t a{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text16.size());
                a += s16_i.rank_limit<2>(pos, symb);
            }
            ankerl::nanobench::doNotOptimizeAway(a);
        });
    }
    #if 1

    SECTION("benchmarking prefix_rank") {
        auto rng = ankerl::nanobench::Rng{};

        auto bench = ankerl::nanobench::Bench{};
        bench.title("prefix_rank")
             .relative(true)
             .batch(100);

        bench.run("prefix_rank-4", [&]() {
            auto symb = rng.bounded(5);
            size_t a{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size());
                a += s4.prefix_rank(pos, symb);
            }
            ankerl::nanobench::doNotOptimizeAway(a);

        });
        bench.run("prefix_rank_limit<2>-4", [&]() {
            auto symb = rng.bounded(5);
            size_t a{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size());
                a += s4.prefix_rank_limit<2>(pos, symb);
            }
            ankerl::nanobench::doNotOptimizeAway(a);

        });
        bench.run("prefix_rank-16", [&]() {
            auto symb = rng.bounded(17);
            size_t a{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text16.size());
                a += s16.prefix_rank(pos, symb);
            }
            ankerl::nanobench::doNotOptimizeAway(a);
        });
        bench.run("prefix_rank_limit<2>-16", [&]() {
            auto symb = rng.bounded(5);
            size_t a{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text16.size());
                a += s16.prefix_rank_limit<2>(pos, symb);
            }
            ankerl::nanobench::doNotOptimizeAway(a);
        });

        bench.run("interleaved prefix_rank-4", [&]() {
            auto symb = rng.bounded(5);
            size_t a{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size());
                a += s4_b.prefix_rank(pos, symb);
            }
            ankerl::nanobench::doNotOptimizeAway(a);

        });
        bench.run("interleaved prefix_rank_limit<2>-4", [&]() {
            auto symb = rng.bounded(5);
            size_t a{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size());
                a += s4_b.prefix_rank_limit(pos, symb);
            }
            ankerl::nanobench::doNotOptimizeAway(a);

        });

        bench.run("interleaved prefix_rank-16", [&]() {
            auto symb = rng.bounded(17);
            size_t a{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text16.size());
                a += s16_b.prefix_rank(pos, symb);
            }
            ankerl::nanobench::doNotOptimizeAway(a);
        });

        bench.run("interleaved prefix_rank_limit<2>-16", [&]() {
            auto symb = rng.bounded(5);
            size_t a{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text16.size());
                a += s16_b.prefix_rank_limit(pos, symb);
            }
            ankerl::nanobench::doNotOptimizeAway(a);
        });
    }
#endif
#if 1
    SECTION("benchmarking prefix_rank_and_rank") {
        auto rng = ankerl::nanobench::Rng{};

        auto bench = ankerl::nanobench::Bench{};
        bench.title("prefix_rank_and_rank")
             .relative(true)
             .batch(100);

        bench.run("prefix_rank_and_rank-4", [&]() {
            auto symb = rng.bounded(4);
            size_t a{}, a2{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size());
                auto [p, p2] = s4.prefix_rank_and_rank(pos, symb);
                a += p;
                a2 += p2;
            }
            ankerl::nanobench::doNotOptimizeAway(a);
            ankerl::nanobench::doNotOptimizeAway(a2);
        });
        bench.run("prefix_rank_and_rank_limit<2>-4", [&]() {
            auto symb = rng.bounded(4);
            size_t a{}, a2{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size());
                auto [p, p2] = s4.prefix_rank_and_rank_limit<2>(pos, symb);
                a += p;
                a2 += p2;

            }
            ankerl::nanobench::doNotOptimizeAway(a);
            ankerl::nanobench::doNotOptimizeAway(a2);
        });
        bench.run("prefix_rank_and_rank-16", [&]() {
            auto symb = rng.bounded(16);
            size_t a{}, a2{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text16.size());
                auto [p, p2] = s16.prefix_rank_and_rank(pos, symb);
                a += p;
                a2 += p2;
            }
            ankerl::nanobench::doNotOptimizeAway(a);
            ankerl::nanobench::doNotOptimizeAway(a2);
        });
        bench.run("prefix_rank_and_rank_limit<2>-16", [&]() {
            auto symb = rng.bounded(4);
            size_t a{}, a2{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text16.size());
                auto [p, p2] = s16.prefix_rank_and_rank_limit<2>(pos, symb);
                a += p;
                a2 += p2;
            }
            ankerl::nanobench::doNotOptimizeAway(a);
            ankerl::nanobench::doNotOptimizeAway(a2);
        });
        bench.run("prefix_rank_and_rank_limit<4>-16", [&]() {
            auto symb = rng.bounded(16);
            size_t a{}, a2{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text16.size());
                auto [p, p2] = s16.prefix_rank_and_rank_limit<4>(pos, symb);
                a += p;
                a2 += p2;
            }
            ankerl::nanobench::doNotOptimizeAway(a);
            ankerl::nanobench::doNotOptimizeAway(a2);
        });

        bench.run("interleaved prefix_rank_and_rank-4", [&]() {
            auto symb = rng.bounded(4);
            size_t a{}, a2{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size());
                auto [p, p2] = s4_b.prefix_rank_and_rank(pos, symb);
                a += p;
                a2 += p2;
            }
            ankerl::nanobench::doNotOptimizeAway(a);
            ankerl::nanobench::doNotOptimizeAway(a2);
        });
        bench.run("interleaved prefix_rank_and_rank_limit<2>-4", [&]() {
            auto symb = rng.bounded(4);
            size_t a{}, a2{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size());
                auto [p, p2] = s4_b.prefix_rank_and_rank_limit(pos, symb);
                a += p;
                a2 += p2;
            }
            ankerl::nanobench::doNotOptimizeAway(a);
            ankerl::nanobench::doNotOptimizeAway(a2);
        });

        bench.run("interleaved prefix_rank_and_rank-16", [&]() {
            auto symb = rng.bounded(16);
            size_t a{}, a2{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text16.size());
                auto [p, p2] = s16_b.prefix_rank_and_rank(pos, symb);
                a += p;
                a2 += p;
            }
            ankerl::nanobench::doNotOptimizeAway(a);
            ankerl::nanobench::doNotOptimizeAway(a2);
        });
        bench.run("interleaved prefix_rank_and_rank_limit<2>-16", [&]() {
            auto symb = rng.bounded(4);
            size_t a{}, a2{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text16.size());
                auto [p, p2] = s16_b.prefix_rank_and_rank_limit(pos, symb);
                a += p;
                a2 += p;
            }
            ankerl::nanobench::doNotOptimizeAway(a);
            ankerl::nanobench::doNotOptimizeAway(a2);
        });
        bench.run("interleaved prefix_rank_and_rank_limit<4>-16", [&]() {
            auto symb = rng.bounded(16);
            size_t a{}, a2{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text16.size());
                auto [p, p2] = s16_b.prefix_rank_and_rank_limit(pos, symb);
                a += p;
                a2 += p;
            }
            ankerl::nanobench::doNotOptimizeAway(a);
            ankerl::nanobench::doNotOptimizeAway(a2);
        });
    }

    SECTION("benchmarking all_ranks") {
        auto rng = ankerl::nanobench::Rng{};

        auto bench = ankerl::nanobench::Bench{};
        bench.title("all_ranks")
             .relative(true)
             .batch(100);

        bench.run("all_ranks-4", [&]() {
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size());
                auto v = s4.all_ranks(pos);
                ankerl::nanobench::doNotOptimizeAway(v);
            }
        });
        bench.run("all_ranks-16", [&]() {
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size());
                auto v = s16.all_ranks(pos);
                ankerl::nanobench::doNotOptimizeAway(v);
            }
        });
        bench.run("interleaved all_ranks-4", [&]() {
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size());
                auto v = s4_b.all_ranks(pos);
                ankerl::nanobench::doNotOptimizeAway(v);
            }
        });
        bench.run("interleaved all_ranks-16", [&]() {
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size());
                auto v = s16_b.all_ranks(pos);
                ankerl::nanobench::doNotOptimizeAway(v);
            }
        });
    }

    SECTION("benchmarking all_ranks_and_prefix_ranks") {
        auto rng = ankerl::nanobench::Rng{};

        auto bench = ankerl::nanobench::Bench{};
        bench.title("all_ranks_and_prefix_ranks")
             .relative(true)
             .batch(100);

        bench.run("all_ranks_and_prefix_ranks-4", [&]() {
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size());
                auto [v, v2] = s4.all_ranks_and_prefix_ranks(pos);
                ankerl::nanobench::doNotOptimizeAway(v);
                ankerl::nanobench::doNotOptimizeAway(v2);
            }
        });
        bench.run("all_ranks_and_prefix_ranks-16", [&]() {
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size());
                auto [v, v2] = s16.all_ranks_and_prefix_ranks(pos);
                ankerl::nanobench::doNotOptimizeAway(v);
                ankerl::nanobench::doNotOptimizeAway(v2);
            }
        });
        bench.run("interleaved all_ranks_and_prefix_ranks-4", [&]() {
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size());
                auto [v, v2] = s4_b.all_ranks_and_prefix_ranks(pos);
                ankerl::nanobench::doNotOptimizeAway(v);
                ankerl::nanobench::doNotOptimizeAway(v2);
            }
        });
        bench.run("interleaved all_ranks_and_prefix_ranks-16", [&]() {
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size());
                auto [v, v2] = s16_b.all_ranks_and_prefix_ranks(pos);
                ankerl::nanobench::doNotOptimizeAway(v);
                ankerl::nanobench::doNotOptimizeAway(v2);
            }
        });
    }


    SECTION("benchmarking all_ranks_dual") {
        auto rng = ankerl::nanobench::Rng{};

        auto bench = ankerl::nanobench::Bench{};
        bench.title("all_ranks_dual")
             .relative(true)
             .batch(100);

        bench.run("all_ranks_dual-4", [&]() {
            size_t a{}, a1{}, a2{}, a3{}, a4{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size()-1);
                auto pos2 = pos + rng.bounded(text4.size()-pos);
                s4.all_ranks_dual(pos, pos2, [&](size_t symb, auto v1, auto v2, auto v3, auto v4) {
                    a += symb;
                    a1 += v1;
                    a2 += v2;
                    a3 += v3;
                    a4 += v4;
                });
            }
            ankerl::nanobench::doNotOptimizeAway(a);
            ankerl::nanobench::doNotOptimizeAway(a1);
            ankerl::nanobench::doNotOptimizeAway(a2);
            ankerl::nanobench::doNotOptimizeAway(a3);
            ankerl::nanobench::doNotOptimizeAway(a4);
        });

        bench.run("all_ranks_dual_limit<2>-4", [&]() {
            size_t a{}, a1{}, a2{}, a3{}, a4{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size()-1);
                auto pos2 = pos + rng.bounded(text4.size()-pos);
                s4.all_ranks_dual_limit<2>(pos, pos2, [&](size_t symb, auto v1, auto v2, auto v3, auto v4) {
                    a += symb;
                    a1 += v1;
                    a2 += v2;
                    a3 += v3;
                    a4 += v4;
                });
            }
            ankerl::nanobench::doNotOptimizeAway(a);
            ankerl::nanobench::doNotOptimizeAway(a1);
            ankerl::nanobench::doNotOptimizeAway(a2);
            ankerl::nanobench::doNotOptimizeAway(a3);
            ankerl::nanobench::doNotOptimizeAway(a4);
        });

        bench.run("all_ranks_dual-16", [&]() {
            size_t a{}, a1{}, a2{}, a3{}, a4{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size()-1);
                auto pos2 = pos + rng.bounded(text4.size()-pos);
                s16.all_ranks_dual(pos, pos2, [&](size_t symb, auto v1, auto v2, auto v3, auto v4) {
                    a += symb;
                    a1 += v1;
                    a2 += v2;
                    a3 += v3;
                    a4 += v4;
                });
            }
            ankerl::nanobench::doNotOptimizeAway(a);
            ankerl::nanobench::doNotOptimizeAway(a1);
            ankerl::nanobench::doNotOptimizeAway(a2);
            ankerl::nanobench::doNotOptimizeAway(a3);
            ankerl::nanobench::doNotOptimizeAway(a4);
        });
        bench.run("all_ranks_dual_limit<2>-16", [&]() {
            size_t a{}, a1{}, a2{}, a3{}, a4{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size()-1);
                auto pos2 = pos + rng.bounded(text4.size()-pos);
                s16.all_ranks_dual_limit<2>(pos, pos2, [&](size_t symb, auto v1, auto v2, auto v3, auto v4) {
                    a += symb;
                    a1 += v1;
                    a2 += v2;
                    a3 += v3;
                    a4 += v4;
                });
            }
            ankerl::nanobench::doNotOptimizeAway(a);
            ankerl::nanobench::doNotOptimizeAway(a1);
            ankerl::nanobench::doNotOptimizeAway(a2);
            ankerl::nanobench::doNotOptimizeAway(a3);
            ankerl::nanobench::doNotOptimizeAway(a4);
        });
        bench.run("all_ranks_dual_limit<4>-16", [&]() {
            size_t a{}, a1{}, a2{}, a3{}, a4{};
            for (size_t i{0}; i < 100; ++i) {
                auto pos = rng.bounded(text4.size()-1);
                auto pos2 = pos + rng.bounded(text4.size()-pos);
                s16.all_ranks_dual_limit<4>(pos, pos2, [&](size_t symb, auto v1, auto v2, auto v3, auto v4) {
                    a += symb;
                    a1 += v1;
                    a2 += v2;
                    a3 += v3;
                    a4 += v4;
                });
            }
            ankerl::nanobench::doNotOptimizeAway(a);
            ankerl::nanobench::doNotOptimizeAway(a1);
            ankerl::nanobench::doNotOptimizeAway(a2);
            ankerl::nanobench::doNotOptimizeAway(a3);
            ankerl::nanobench::doNotOptimizeAway(a4);
        });
    }
    #endif
    #endif
}
