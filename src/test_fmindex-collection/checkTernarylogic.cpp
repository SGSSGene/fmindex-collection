//SPDX-FileCopyrightText: 2024 Simon Gene Gottlieb
//SPDX-License-Identifier: BSD-3-Clause

#include <bit>
#include <catch2/catch_all.hpp>
#include <nanobench.h>
#include <fmindex-collection/ternarylogic.h>


TEST_CASE("testing ternary Logic functions are equivalent", "[ternary_logic]") {
    auto check_equality = [&]<size_t Bits>() {
        INFO(Bits);
        auto a = std::bitset<Bits>{0b00001111};
        auto b = std::bitset<Bits>{0b00110011};
        auto c = std::bitset<Bits>{0b01010101};

        auto dumbTernary = [&]<size_t N>(size_t x, std::bitset<N> const& a, std::bitset<N> const& b, std::bitset<N> const& c) {
            auto dumbBinaryTernary = [&](size_t x, bool a, bool b, bool c) -> bool {
                uint64_t v = (a<<2) + (b<<1) + c;
                return x & (1<<v);
            };

            auto r = std::bitset<N>{};
            for (size_t i{0}; i < N; ++i) {
                r[i] = dumbBinaryTernary(x, a[i], b[i], c[i]);
            }
            return r;
        };

        fmc::for_constexpr<0, 256>([&]<size_t I>() {
            INFO(I);
            auto v0 = dumbTernary(I, a, b, c);
            auto v1 = fmc::ternarylogic_v1<I, Bits>(a, b, c);
            auto v2 = fmc::ternarylogic_v2<Bits>(I, a, b, c);
            auto v3 = fmc::ternarylogic_v3<Bits>(I, a, b, c);
//            auto v4 = and_and<I>(a, b, c);
//            auto v5 = and_and_run(I, a, b, c);
//            auto v6 = fast_ternarylogic(I, a, b, c);

            CHECK(v0 == v1);
            CHECK(v0 == v2);
            CHECK(v0 == v3);
//            CHECK(v0 == v4);
//            CHECK(v0 == v5);
//            CHECK(v0 == v6);
        });
    };

    SECTION("check simde and and_and equal") {
        check_equality.operator()<64>();
        check_equality.operator()<128>();
        check_equality.operator()<256>();
        check_equality.operator()<512>();
        check_equality.operator()<1024>();
        check_equality.operator()<2048>();
        check_equality.operator()<4096>();
        check_equality.operator()<8192>();
        check_equality.operator()<16384>();
        check_equality.operator()<32768>();
        check_equality.operator()<65536>();
    }
}

TEST_CASE("mark_exact - test to ternary logic function", "[ternary_logic][mark_exact]") {

    auto check_equality = [&]<size_t Bits>() {
        auto a = std::bitset<Bits>{0b00001111};
        auto b = std::bitset<Bits>{0b00110011};
        auto c = std::bitset<Bits>{0b01010101};

        auto f = [&a, &b, &c]<size_t N1, size_t N2>() {
            auto expected = fmc::ternarylogic_v3<Bits>(N2, a, b, c);
            auto r2       = fmc::mark_exact_v2<Bits>(N1, a, b, c);
            auto r3       = fmc::mark_exact_v3<Bits>(N1, a, b, c);
            auto r4       = fmc::mark_exact_v4<Bits>(N1, a, b, c);
            auto rfast    = fmc::mark_exact_fast<Bits>(N1, a, b, c);
            auto r_all    = fmc::mark_exact_all<Bits>(a, b, c)[N1];
            CHECK(expected == r2);
            CHECK(expected == r3);
            CHECK(expected == r4);
            CHECK(expected == rfast);
            CHECK(expected == r_all);
        };
        f.template operator()<0,   1>();
        f.template operator()<1,   2>();
        f.template operator()<2,   4>();
        f.template operator()<3,   8>();
        f.template operator()<4,  16>();
        f.template operator()<5,  32>();
        f.template operator()<6,  64>();
        f.template operator()<7, 128>();
    };
    SECTION("check ternary and mark_exact are equal") {
        check_equality.operator()<64>();
        check_equality.operator()<128>();
        check_equality.operator()<256>();
        check_equality.operator()<512>();
        check_equality.operator()<1024>();
        check_equality.operator()<2048>();
        check_equality.operator()<4096>();
        check_equality.operator()<8192>();
        check_equality.operator()<16384>();
        check_equality.operator()<32768>();
        check_equality.operator()<65536>();
    }
}

TEST_CASE("mark_exact - testing and benchmarking", "[ternary_logic][!benchmark][mark_exact]") {
    auto rng = ankerl::nanobench::Rng{};
    SECTION("test all of it") {
        auto benchmark = [&]<size_t Bits>() {
            using BS = std::bitset<Bits>;
            auto bitsets = std::vector<BS>{};
            #ifdef NDEBUG
            bitsets.reserve(10'000'000*64 / Bits);
            #else
            bitsets.reserve(1'000*64 / Bits);
            if (bitsets.capacity() == 0) {
                bitsets.reserve(1);
            }
            #endif
            for (size_t i{0}; i < bitsets.capacity(); ++i) {
                auto b = BS{};
                for (size_t j{0}; j < Bits/64; ++j) {
                    b = (b << 64) | BS{rng()};
                }
                bitsets.push_back(b);
            }

            auto bits_as_str = std::to_string(Bits);

            auto bench = ankerl::nanobench::Bench{};
            bench.title("mark exact " + bits_as_str)
                 .relative(true)
                 .epochs(100)
                 .batch(256);

            auto i1 = rng.bounded(bitsets.size());
            auto i2 = rng.bounded(bitsets.size());
            auto i3 = rng.bounded(bitsets.size());
            auto offset = rng.bounded(8);

            bench.run("mark_exact_v2 " + bits_as_str, [&]() {
                for (size_t i{0}; i < 32; ++i) {
                fmc::for_constexpr<0, 8>([&]<size_t I>() {
                    auto v = fmc::mark_exact_v2((I+offset)%8, bitsets[i1], bitsets[i2], bitsets[i3]);
                    ankerl::nanobench::doNotOptimizeAway(v);
                });
                }
            });
            bench.run("mark_exact_v3 " + bits_as_str, [&]() {
                for (size_t i{0}; i < 32; ++i) {
                fmc::for_constexpr<0, 8>([&]<size_t I>() {
                    auto v = fmc::mark_exact_v3((I+offset)%8, bitsets[i1], bitsets[i2], bitsets[i3]);
                    ankerl::nanobench::doNotOptimizeAway(v);
                });
                }
            });
            bench.run("mark_exact_v4 " + bits_as_str, [&]() {
                for (size_t i{0}; i < 32; ++i) {
                fmc::for_constexpr<0, 8>([&]<size_t I>() {
                    auto v = fmc::mark_exact_v4((I+offset)%8, bitsets[i1], bitsets[i2], bitsets[i3]);
                    ankerl::nanobench::doNotOptimizeAway(v);
                });
                }
            });

            bench.run("mark_exact_fast " + bits_as_str, [&]() {
                for (size_t i{0}; i < 32; ++i) {
                fmc::for_constexpr<0, 8>([&]<size_t I>() {
                    auto v = fmc::mark_exact_fast((I+offset)%8, bitsets[i1], bitsets[i2], bitsets[i3]);
                    ankerl::nanobench::doNotOptimizeAway(v);
                });
                }
            });
            bench.run("mark_exact_all " + bits_as_str, [&]() {
                for (size_t i{0}; i < 32; ++i) {
                auto v = fmc::mark_exact_all(bitsets[i1], bitsets[i2], bitsets[i3]);
                fmc::for_constexpr<0, 8>([&]<size_t I>() {
                    ankerl::nanobench::doNotOptimizeAway(v[(I+offset)%8]);
                });
                }
            });

        };
        benchmark.operator()<32>();
        benchmark.operator()<64>();
        benchmark.operator()<128>();
        benchmark.operator()<256>();
        benchmark.operator()<512>();
        benchmark.operator()<1024>();
        benchmark.operator()<2048>();
        benchmark.operator()<4096>();
        benchmark.operator()<8192>();
        benchmark.operator()<16384>();
        benchmark.operator()<32768>();
        benchmark.operator()<65536>();
    }
}

TEST_CASE("mark_exact_or_less - test to ternary logic function", "[ternary_logic][mark_exact_or_less]") {
    auto check_equality = [&]<size_t Bits>() {
        auto a = std::bitset<Bits>{0b00001111};
        auto b = std::bitset<Bits>{0b00110011};
        auto c = std::bitset<Bits>{0b01010101};

        auto f = [&a, &b, &c]<size_t N1, size_t N2>() {
            auto expected = fmc::ternarylogic_v3<Bits>(N2, a, b, c);
            auto r2       = fmc::mark_exact_or_less_v2<Bits>(N1, a, b, c);
            auto r3       = fmc::mark_exact_or_less_v3<Bits>(N1, a, b, c);
            auto r4       = fmc::mark_exact_or_less_v4<Bits>(N1, a, b, c);
            auto rfast    = fmc::mark_exact_or_less_fast<Bits>(N1, a, b, c);
            auto r_all    = fmc::mark_exact_or_less_all<Bits>(a, b, c)[N1];
            CHECK(expected == r2);
            CHECK(expected == r3);
            CHECK(expected == r4);
            CHECK(expected == rfast);
            CHECK(expected == r_all);
        };
        f.template operator()<0,   1>();
        f.template operator()<1,   3>();
        f.template operator()<2,   7>();
        f.template operator()<3,  15>();
        f.template operator()<4,  31>();
        f.template operator()<5,  63>();
        f.template operator()<6, 127>();
        f.template operator()<7, 255>();

    };
    SECTION("check mark_exact_or_less and ternary are equal") {
        check_equality.operator()<64>();
        check_equality.operator()<128>();
        check_equality.operator()<256>();
        check_equality.operator()<512>();
        check_equality.operator()<1024>();
        check_equality.operator()<2048>();
        check_equality.operator()<4096>();
        check_equality.operator()<8192>();
        check_equality.operator()<16384>();
        check_equality.operator()<32768>();
        check_equality.operator()<65536>();
    }
}

//!WORKAROUND !TODO not running on msvc because 14.44 is crashing with internal compiler error (14.43 is working fine)
 #if !defined(_MSC_VER) || defined(__clang__) // not MSVC
TEST_CASE("mark_exact_or_less_large - test to ternary logic function", "[ternary_logic][mark_exact_or_less_large]") {
    auto check_equality = [&]<size_t Bits, size_t BitWidth=4>() {
        for (size_t loop{0}; loop < 1<<BitWidth; loop += Bits) {
            auto b = std::array<std::bitset<Bits>, BitWidth>{};
            for (size_t i{loop}; i < loop + Bits; ++i) {
                for (size_t j{0}; j < b.size(); ++j) {
                    b[j].set(i-loop, (i>>j)&1);
                }
            }

            auto f = [&]<size_t CheckValue>() {
                auto dumbTernary = [&]() {
                    auto reconstructValue = [&](size_t idx) {
                        size_t v{};
                        for (size_t i{0}; i<b.size(); ++i) {
                            v = v | (b[i].test(idx) << i);
                        }
                        return v;
                    };
                    auto r = std::bitset<Bits>{};
                    for (size_t i{0}; i < Bits; ++i) {
                        auto v = reconstructValue(i);
                        r.set(i, v <= CheckValue);
                    }
                    return r;
                };
                INFO(b[0]);
                INFO(CheckValue);

                auto expected = dumbTernary();
                auto r        = fmc::mark_exact_or_less_large<Bits>(CheckValue, b);
                auto r_exact  = [&]() {
                    auto r = fmc::mark_exact_large(0, b);
                    for (size_t i{1}; i <= CheckValue; ++i) {
                        r = r | fmc::mark_exact_large(i, b);
                    }
                    return r;
                }();
                CHECK(expected == r);
                CHECK(expected == r_exact);

                if constexpr (BitWidth == 3) {
                    auto r_ternary = fmc::mark_exact_or_less_v3<Bits>(CheckValue, b[2], b[1], b[0]);
                    CHECK(expected == r_ternary);
                }
            };
            fmc::for_constexpr<0, (1<<BitWidth)>([&]<size_t I>() {
                f.template operator()<I>();
            });
        }
    };
    SECTION("check mark_exact_or_less_large and ternary are equal") {
        fmc::for_constexpr<1, 9>([&]<size_t I>() {
            INFO(I << "info I: ");
            check_equality.operator()<64, I>();
            check_equality.operator()<256, I>();
            check_equality.operator()<512>();
            check_equality.operator()<1024>();
            check_equality.operator()<2048>();
            check_equality.operator()<4096>();
            check_equality.operator()<8192>();
            check_equality.operator()<16384>();
            check_equality.operator()<32768>();
            check_equality.operator()<65536>();
        });
    }
    SECTION("check specific errornous configuration") {
        uint8_t symb = 32;
        auto arr = std::array<std::bitset<64>, 8>{};
        for (auto p : {1, 4, 6, 7}) arr[0].set(p);
        for (auto p : {4, 6}) arr[1].set(p);
        for (auto p : {2, 3, 4, 6, 7, 8, 9}) arr[2].set(p);
        for (auto p : {0, 2, 3, 4, 8}) arr[3].set(p);
        for (auto p : {6, 9}) arr[4].set(p);
        for (auto p : {1, 2, 3, 4, 5, 7, 8, 9}) arr[5].set(p);
        for (auto p : {0, 1, 2, 3, 4, 6, 7, 8, 9}) arr[6].set(p);
//        for (auto p : {}) arr[7].set(p);

        auto v = fmc::mark_exact_or_less_large(symb, arr);
        auto mask = fmc::mark_exact_large(0, arr);
        for (uint64_t i{1}; i <= symb; ++i) {
            mask |= fmc::mark_exact_large(i, arr);
        }
        CHECK(v == mask);
    }
}
#endif
#if 1

TEST_CASE("mark_exact_or_less - testing and benchmarking", "[ternary_logic][!benchmark][mark_exact_or_less]") {
    auto rng = ankerl::nanobench::Rng{};
    SECTION("test all of it") {
        auto benchmark = [&]<size_t Bits>() {
            using BS = std::bitset<Bits>;
            auto bitsets = std::vector<BS>{};
            #ifdef NDEBUG
            bitsets.reserve(10'000'000*64 / Bits);
            #else
            bitsets.reserve(1'000*64 / Bits);
            if (bitsets.capacity() == 0) {
                bitsets.reserve(1);
            }
            #endif
            for (size_t i{0}; i < bitsets.capacity(); ++i) {
                auto b = BS{};
                for (size_t j{0}; j < Bits/64; ++j) {
                    b = (b << 64) | BS{rng()};
                }
                bitsets.push_back(b);
            }

            auto bits_as_str = std::to_string(Bits);

            auto bench = ankerl::nanobench::Bench{};
            bench.title("mark exact " + bits_as_str)
                 .relative(true)
                 .epochs(100)
                 .batch(256);

            auto i1 = rng.bounded(bitsets.size());
            auto i2 = rng.bounded(bitsets.size());
            auto i3 = rng.bounded(bitsets.size());
            auto offset = rng.bounded(8);

            bench.run("mark_exact_or_less_v2 " + bits_as_str, [&]() {
                for (size_t i{0}; i < 32; ++i) {
                fmc::for_constexpr<0, 8>([&]<size_t I>() {
                    auto v = fmc::mark_exact_or_less_v2((I+offset)%8, bitsets[i1], bitsets[i2], bitsets[i3]);
                    ankerl::nanobench::doNotOptimizeAway(v);
                });
                }
            });
            bench.run("mark_exact_or_less_v3 " + bits_as_str, [&]() {
                for (size_t i{0}; i < 32; ++i) {
                fmc::for_constexpr<0, 8>([&]<size_t I>() {
                    auto v = fmc::mark_exact_or_less_v3((I+offset)%8, bitsets[i1], bitsets[i2], bitsets[i3]);
                    ankerl::nanobench::doNotOptimizeAway(v);
                });
                }
            });
            bench.run("mark_exact_or_less_v4 " + bits_as_str, [&]() {
                for (size_t i{0}; i < 32; ++i) {
                fmc::for_constexpr<0, 8>([&]<size_t I>() {
                    auto v = fmc::mark_exact_or_less_v4((I+offset)%8, bitsets[i1], bitsets[i2], bitsets[i3]);
                    ankerl::nanobench::doNotOptimizeAway(v);
                });
                }
            });

            bench.run("mark_exact_or_less_fast " + bits_as_str, [&]() {
                for (size_t i{0}; i < 32; ++i) {
                fmc::for_constexpr<0, 8>([&]<size_t I>() {
                    auto v = fmc::mark_exact_or_less_fast((I+offset)%8, bitsets[i1], bitsets[i2], bitsets[i3]);
                    ankerl::nanobench::doNotOptimizeAway(v);
                });
                }
            });

            bench.run("mark_exact_or_less_all " + bits_as_str, [&]() {
                for (size_t i{0}; i < 32; ++i) {
                    auto v = fmc::mark_exact_or_less_all(bitsets[i1], bitsets[i2], bitsets[i3]);
                    fmc::for_constexpr<0, 8>([&]<size_t I>() {
                        ankerl::nanobench::doNotOptimizeAway(v[(I+offset)%8]);
                    });
                }
            });

        };
        benchmark.operator()<32>();
        benchmark.operator()<64>();
        benchmark.operator()<128>();
        benchmark.operator()<256>();
        benchmark.operator()<512>();
        benchmark.operator()<1024>();
        benchmark.operator()<2048>();
        benchmark.operator()<4096>();
        benchmark.operator()<8192>();
        benchmark.operator()<16384>();
        benchmark.operator()<32768>();
        benchmark.operator()<65536>();
    }
}
#endif
