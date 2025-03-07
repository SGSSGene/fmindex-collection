//SPDX-FileCopyrightText: 2024 Simon Gene Gottlieb
//SPDX-License-Identifier: BSD-3-Clause

#include <bit>
#include <catch2/catch_all.hpp>
#include <nanobench.h>
#include <fmindex-collection/bitset_popcount.h>
#include <fmindex-collection/ternarylogic.h>


/*
#include <simde/x86/avx512.h>
#include <simde/x86/sse2.h>

template <size_t T, size_t O>
auto convert(std::bitset<O> const& _o) -> std::bitset<T> {
    using A1 = std::array<uint64_t, (O+63)/64>;
    auto const& s = reinterpret_cast<A1 const&>(_o);

    auto t = std::array<uint64_t, (T+63)/64>{};
    for (size_t i{0}; i < std::min(s.size(), t.size()); ++i) {
        t[i] = s[i];
    }
    return std::bit_cast<std::bitset<T>>(t);
}


template <size_t Shift, size_t N>
auto and_and(std::bitset<N> const& _a, std::bitset<N> const& _b, std::bitset<N> const& _c) -> std::bitset<N> {
    if constexpr (N == 1024 || N == 2048) {
        (void)_a;
        (void)_b;
        (void)_c;
        return {};
    } else if constexpr (N == 512) {
        auto a = simde_mm512_loadu_epi64(&_a);
        auto b = simde_mm512_loadu_epi64(&_b);
        auto c = simde_mm512_loadu_epi64(&_c);
        auto r = simde_mm512_ternarylogic_epi64(a, b, c, Shift);
        auto _r = std::bitset<N>{};
        simde_mm512_storeu_epi64(&_r, r);
        return _r;
    } else if constexpr (N == 256) {
        auto a = simde_mm256_loadu_epi64(&_a);
        auto b = simde_mm256_loadu_epi64(&_b);
        auto c = simde_mm256_loadu_epi64(&_c);
        auto r = simde_mm256_ternarylogic_epi64(a, b, c, Shift);
        auto _r = std::bitset<N>{};
        simde_mm256_storeu_epi64(&_r, r);
        return _r;
    } else if constexpr (N == 128) {
        auto r = and_and<Shift>(convert<256>(_a), convert<256>(_b), convert<256>(_c));
        return convert<128>(r);
    } else if constexpr (N == 64) {
        auto r = and_and<Shift>(convert<256>(_a), convert<256>(_b), convert<256>(_c));
        return convert<64>(r);
    } else {
        []<bool b = false> {
            static_assert(b, "Only values of 64, 128, 256 and 512 supported");
        };
    }
}

template <size_t N>
auto and_and_run(size_t shift, std::bitset<N> const& _a, std::bitset<N> const& _b, std::bitset<N> const& _c) -> std::bitset<N> {
    if constexpr (N == 1024 || N == 2048) {
        (void)shift;
        (void)_a;
        (void)_b;
        (void)_c;
        return {};
    } else if constexpr (N == 512) {
        auto a = simde_mm512_loadu_epi64(&_a);
        auto b = simde_mm512_loadu_epi64(&_b);
        auto c = simde_mm512_loadu_epi64(&_c);
        auto r = simde_mm512_ternarylogic_epi64(a, b, c, shift);
        auto _r = std::bitset<N>{};
        simde_mm512_storeu_epi64(&_r, r);
        return _r;
    } else if constexpr (N == 256) {
        auto a = simde_mm256_loadu_epi64(&_a);
        auto b = simde_mm256_loadu_epi64(&_b);
        auto c = simde_mm256_loadu_epi64(&_c);
        auto r = simde_mm256_ternarylogic_epi64(a, b, c, shift);
        auto _r = std::bitset<N>{};
        simde_mm256_storeu_epi64(&_r, r);
        return _r;
    } else if constexpr (N == 128) {
        auto r = and_and_run(shift, convert<256>(_a), convert<256>(_b), convert<256>(_c));
        return convert<128>(r);
    } else if constexpr (N == 64) {
        auto r = and_and_run(shift, convert<256>(_a), convert<256>(_b), convert<256>(_c));
        return convert<64>(r);
    } else {
        []<bool b = false> {
            static_assert(b, "Only values of 64, 128, 256 and 512 supported");
        };
    }
}

template <size_t N>
auto fast_ternarylogic(size_t shift, std::bitset<N> const& _a, std::bitset<N> const& _b, std::bitset<N> const& _c) -> std::bitset<N> {
    if constexpr (N >= 1024) {
        return ternarylogic_v3<N>(shift, _a, _b, _c);
    } else if constexpr (N == 512) {
        auto a = simde_mm512_loadu_epi64(&_a);
        auto b = simde_mm512_loadu_epi64(&_b);
        auto c = simde_mm512_loadu_epi64(&_c);
        auto r = simde_mm512_ternarylogic_epi64(a, b, c, shift);
        auto _r = std::bitset<N>{};
        simde_mm512_storeu_epi64(&_r, r);
        return _r;
    } else if constexpr (N == 256) {
        auto a = simde_mm256_loadu_epi64(&_a);
        auto b = simde_mm256_loadu_epi64(&_b);
        auto c = simde_mm256_loadu_epi64(&_c);
        auto r = simde_mm256_ternarylogic_epi64(a, b, c, shift);
        auto _r = std::bitset<N>{};
        simde_mm256_storeu_epi64(&_r, r);
        return _r;
    } else if constexpr (N == 128) {
        return ternarylogic_v3<N>(shift, _a, _b, _c);
    } else if constexpr (N == 64) {
        return ternarylogic_v3<N>(shift, _a, _b, _c);
    } else {
        []<bool b = false> {
            static_assert(b, "Only values of 64, 128, 256 and 512 supported");
        };
    }
}*/

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

        fmindex_collection::for_constexpr<0, 256>([&]<size_t I>() {
            INFO(I);
            auto v0 = dumbTernary(I, a, b, c);
            auto v1 = fmindex_collection::ternarylogic_v1<I, Bits>(a, b, c);
            auto v2 = fmindex_collection::ternarylogic_v2<Bits>(I, a, b, c);
            auto v3 = fmindex_collection::ternarylogic_v3<Bits>(I, a, b, c);
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

TEST_CASE("mark_exact - test to ternary logic function", "[mark_exact]") {

    auto check_equality = [&]<size_t Bits>() {
        auto a = std::bitset<Bits>{0b00001111};
        auto b = std::bitset<Bits>{0b00110011};
        auto c = std::bitset<Bits>{0b01010101};

        auto f = [&a, &b, &c]<size_t N1, size_t N2>() {
            auto expected = fmindex_collection::ternarylogic_v3<Bits>(N2, a, b, c);
            auto r2       = fmindex_collection::mark_exact_v2<Bits>(N1, a, b, c);
            auto r3       = fmindex_collection::mark_exact_v3<Bits>(N1, a, b, c);
            auto r4       = fmindex_collection::mark_exact_v4<Bits>(N1, a, b, c);
            auto rfast    = fmindex_collection::mark_exact_fast<Bits>(N1, a, b, c);
            CHECK(expected == r2);
            CHECK(expected == r3);
            CHECK(expected == r4);
            CHECK(expected == rfast);
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

TEST_CASE("mark_exact - testing and benchmarking", "[mark_exact][!benchmark]") {
    auto rng = ankerl::nanobench::Rng{};
    SECTION("test all of it") {
        auto benchmark = [&]<size_t Bits>() {
            using BS = std::bitset<Bits>;
            auto bitsets = std::vector<BS>{};
            bitsets.reserve(10'000'000*64 / Bits);
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
                fmindex_collection::for_constexpr<0, 8>([&]<size_t I>() {
                    auto v = fmindex_collection::mark_exact_v2((I+offset)%8, bitsets[i1], bitsets[i2], bitsets[i3]);
                    ankerl::nanobench::doNotOptimizeAway(v);
                });
                }
            });
            bench.run("mark_exact_v3 " + bits_as_str, [&]() {
                for (size_t i{0}; i < 32; ++i) {
                fmindex_collection::for_constexpr<0, 8>([&]<size_t I>() {
                    auto v = fmindex_collection::mark_exact_v3((I+offset)%8, bitsets[i1], bitsets[i2], bitsets[i3]);
                    ankerl::nanobench::doNotOptimizeAway(v);
                });
                }
            });
            bench.run("mark_exact_v4 " + bits_as_str, [&]() {
                for (size_t i{0}; i < 32; ++i) {
                fmindex_collection::for_constexpr<0, 8>([&]<size_t I>() {
                    auto v = fmindex_collection::mark_exact_v4((I+offset)%8, bitsets[i1], bitsets[i2], bitsets[i3]);
                    ankerl::nanobench::doNotOptimizeAway(v);
                });
                }
            });

            bench.run("mark_exact_fast " + bits_as_str, [&]() {
                for (size_t i{0}; i < 32; ++i) {
                fmindex_collection::for_constexpr<0, 8>([&]<size_t I>() {
                    auto v = fmindex_collection::mark_exact_fast((I+offset)%8, bitsets[i1], bitsets[i2], bitsets[i3]);
                    ankerl::nanobench::doNotOptimizeAway(v);
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

TEST_CASE("mark_exact_or_less - test to ternary logic function", "[mark_exact_or_less]") {
    auto check_equality = [&]<size_t Bits>() {
        auto a = std::bitset<Bits>{0b00001111};
        auto b = std::bitset<Bits>{0b00110011};
        auto c = std::bitset<Bits>{0b01010101};

        auto f = [&a, &b, &c]<size_t N1, size_t N2>() {
            auto expected = fmindex_collection::ternarylogic_v3<Bits>(N2, a, b, c);
            auto r2       = fmindex_collection::mark_exact_or_less_v2<Bits>(N1, a, b, c);
            auto r3       = fmindex_collection::mark_exact_or_less_v3<Bits>(N1, a, b, c);
            auto r4       = fmindex_collection::mark_exact_or_less_v4<Bits>(N1, a, b, c);
            auto rfast    = fmindex_collection::mark_exact_or_less_fast<Bits>(N1, a, b, c);
            CHECK(expected == r2);
            CHECK(expected == r3);
            CHECK(expected == r4);
            CHECK(expected == rfast);
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

TEST_CASE("mark_exact_or_less_or_less - testing and benchmarking", "[mark_exact_or_less][!benchmark]") {
    auto rng = ankerl::nanobench::Rng{};
    SECTION("test all of it") {
        auto benchmark = [&]<size_t Bits>() {
            using BS = std::bitset<Bits>;
            auto bitsets = std::vector<BS>{};
            bitsets.reserve(10'000'000*64 / Bits);
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
                fmindex_collection::for_constexpr<0, 8>([&]<size_t I>() {
                    auto v = fmindex_collection::mark_exact_or_less_v2((I+offset)%8, bitsets[i1], bitsets[i2], bitsets[i3]);
                    ankerl::nanobench::doNotOptimizeAway(v);
                });
                }
            });
            bench.run("mark_exact_or_less_v3 " + bits_as_str, [&]() {
                for (size_t i{0}; i < 32; ++i) {
                fmindex_collection::for_constexpr<0, 8>([&]<size_t I>() {
                    auto v = fmindex_collection::mark_exact_or_less_v3((I+offset)%8, bitsets[i1], bitsets[i2], bitsets[i3]);
                    ankerl::nanobench::doNotOptimizeAway(v);
                });
                }
            });
            bench.run("mark_exact_or_less_v4 " + bits_as_str, [&]() {
                for (size_t i{0}; i < 32; ++i) {
                fmindex_collection::for_constexpr<0, 8>([&]<size_t I>() {
                    auto v = fmindex_collection::mark_exact_or_less_v4((I+offset)%8, bitsets[i1], bitsets[i2], bitsets[i3]);
                    ankerl::nanobench::doNotOptimizeAway(v);
                });
                }
            });

            bench.run("mark_exact_or_less_fast " + bits_as_str, [&]() {
                for (size_t i{0}; i < 32; ++i) {
                fmindex_collection::for_constexpr<0, 8>([&]<size_t I>() {
                    auto v = fmindex_collection::mark_exact_or_less_fast((I+offset)%8, bitsets[i1], bitsets[i2], bitsets[i3]);
                    ankerl::nanobench::doNotOptimizeAway(v);
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
