// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include <catch2/catch_all.hpp>
#include <fmindex-collection/bitset_popcount.h>

#include <fstream>
#include <nanobench.h>

TEST_CASE("check if signed_rshift_and_count works", "[signed_rshift_and_count]") {
    SECTION("all ones") {
        auto b = std::bitset<64>{};
        b.flip();
        for (size_t i{0}; i < 65; ++i) {
            CHECK(fmc::signed_rshift_and_count(b, i) == 64-i);
        }

        for (size_t i{64}; i < 128; ++i) {
            CHECK(fmc::signed_rshift_and_count(b, i) == i-64);
        }
    }
}
TEST_CASE("count_first_or_last_n_bits", "[count_bits]") {
    SECTION("all ones") {
        auto b = std::bitset<64>{};
        b.flip();
        for (size_t i{0}; i < 65; ++i) {
            CHECK(fmc::skip_first_or_last_n_bits_and_count(b, i) == 64-i);
        }

        for (size_t i{64}; i <= 128; ++i) {
            CHECK(fmc::skip_first_or_last_n_bits_and_count(b, i) == i-64);
        }
    }
    SECTION("single one on the first bit") {
        auto b = std::bitset<64>{};
        b[0] = true;
        CHECK(fmc::skip_first_or_last_n_bits_and_count(b, 0) == 1);
        for (size_t i{1}; i <= 64; ++i) {
            CHECK(fmc::skip_first_or_last_n_bits_and_count(b, i) == 0);
        }
        for (size_t i{65}; i <= 128; ++i) {
            CHECK(fmc::skip_first_or_last_n_bits_and_count(b, i) == 1);
        }
    }

    SECTION("single one on the last bit") {
        auto b = std::bitset<64>{};
        b[63] = true;
        for (size_t i{0}; i < 64; ++i) {
            CHECK(fmc::skip_first_or_last_n_bits_and_count(b, i) == 1);
        }
        for (size_t i{64}; i < 128; ++i) {
            CHECK(fmc::skip_first_or_last_n_bits_and_count(b, i) == 0);
        }
        CHECK(fmc::skip_first_or_last_n_bits_and_count(b, 128) == 1);
    }
}

TEST_CASE("benchmark skipe_first_or_last_n_bits", "[misc][!benchmark]") {
    auto rng = ankerl::nanobench::Rng{};

    SECTION("test shift and count for uint64_t") {
        uint64_t bitset{};
        for (size_t i{0}; i < 64; ++i) {
            bitset = (bitset<<1) | rng.bounded(2);
        }
        auto bench = ankerl::nanobench::Bench{};
        bench.minEpochIterations(10'000).epochs(100);
        bench.relative(true);
            bench.run("bit shift and count - uint64_t" , [&]() {
            auto s = rng.bounded(64);
            auto v = std::popcount(bitset << s);
            ankerl::nanobench::doNotOptimizeAway(v);
        });
    }

    SECTION("test shift and count") {
        auto test = [&]<size_t Bits>() {

            auto bitset = std::bitset<Bits>{};
            for (size_t i{0}; i < Bits; ++i) {
                bitset[i] = rng.bounded(2);
            }
            auto bench = ankerl::nanobench::Bench{};
            bench.minEpochIterations(10'000).epochs(100);
            bench.relative(true);

            auto bits_str = std::to_string(Bits);
            bench.run("bit shift and count " + bits_str, [&]() {
                auto s = rng.bounded(Bits+1);
/*                for (size_t i{0}; i < Bits; ++i) {
                    auto v = (bitset << (s+i)%Bits).count();
                    ankerl::nanobench::doNotOptimizeAway(v);
                }*/
                auto v = (bitset << s).count();
                ankerl::nanobench::doNotOptimizeAway(v);

            });

            bench.run("bit mask and count " + bits_str, [&]() {
                auto s = rng.bounded(Bits+1);
/*                for (size_t i{0}; i < Bits; ++i) {
                    auto v = fmc::skip_first_or_last_n_bits_and_count(bitset, (s+i)%Bits);
                    ankerl::nanobench::doNotOptimizeAway(v);
                }*/
                auto v = fmc::skip_first_or_last_n_bits_and_count(bitset, s);
                ankerl::nanobench::doNotOptimizeAway(v);

            });
        };
        test.template operator()<64>();
        test.template operator()<128>();
        test.template operator()<256>();
        test.template operator()<512>();
        test.template operator()<1024>();
        test.template operator()<2048>();
    }

 /*   SECTION("and and and") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("and and and")
             .relative(true);

        bench.run("naive and 512", [&]() {
            auto i1 = rng.bounded(bitsets.size());
            auto i2 = rng.bounded(bitsets.size());
            auto i3 = rng.bounded(bitsets.size());
            auto v = bitsets[i1] & bitsets[i2] & bitsets[i3];
            ankerl::nanobench::doNotOptimizeAway(v);
        });

        bench.run("naive and_and 512", [&]() {
            auto i1 = rng.bounded(bitsets.size());
            auto i2 = rng.bounded(bitsets.size());
            auto i3 = rng.bounded(bitsets.size());
            auto v = and_and<128>(bitsets[i1], bitsets[i2], bitsets[i3]);
            ankerl::nanobench::doNotOptimizeAway(v);
        });
    }*/
}
