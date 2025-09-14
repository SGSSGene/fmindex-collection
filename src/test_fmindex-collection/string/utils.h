// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#pragma once

#include <catch2/catch_all.hpp>
#include <fmindex-collection/utils.h>
#include <fstream>
#include <nanobench.h>

#include "../BenchSize.h"

#if __has_include(<cxxabi.h>)
#include <cxxabi.h>
namespace {
template <typename T>
auto getName() {
    int     status;
    auto realname = abi::__cxa_demangle(typeid(T).name(), nullptr, nullptr, &status);
    auto str = std::string{realname};
    std::free(realname);
    return str;
}
}
#else
#include <reflect>
namespace {
template <typename T>
auto getName() {
    return std::string{reflect::type_name<T>()};
}
}
#endif

namespace {

template <typename... T>
void call_with_templates(auto f, std::variant<T...> const&) {
    auto g = [&]<typename T2>() {
        if constexpr (!std::same_as<T2, std::monostate>) {
            f.template operator()<T2>();
        }
    };
    (g.template operator()<T>(), ...);
}

template <typename T>
void call_with_templates(auto f) {
    call_with_templates(f, T{});
}

template <template <size_t> class... T>
struct Variant {};

template <size_t>
struct Delimiter {};

template <template <size_t> class... T>
void call_with_templates(auto f, Variant<T...> const&) {
    auto g = [&]<template <size_t> typename T2>() {
        if constexpr (!std::same_as<T2<4>, Delimiter<4>>) {
            f.template operator()<T2>();
        }
    };
    (g.template operator()<T>(), ...);
}

template <typename ...Ts>
struct AppendImpl;

template <typename... T1, typename... Ts>
struct AppendImpl<std::variant<T1...>, Ts...> {
    using type = std::variant<T1..., Ts...>;
};


template <typename T1, typename ...Ts>
using Append = AppendImpl<T1, Ts...>::type;

template <template <size_t, size_t, size_t> class String, size_t l1, size_t l0>
struct Instance {
    template <size_t TSigma>
    using Type = String<TSigma, l1, l0>;
};



template <size_t min, size_t range>
auto generateText(size_t length) -> std::vector<uint8_t> {
    auto rng = ankerl::nanobench::Rng{};

    auto text = std::vector<uint8_t>{};
    for (size_t i{0}; i<length; ++i) {
        text.push_back(rng.bounded(range) + min);
    }
    return text;
}

template <size_t min, size_t range>
auto generateText() -> std::vector<uint8_t> const& {
    static auto text = []() -> std::vector<uint8_t> {
        auto rng = ankerl::nanobench::Rng{};

        // generates string with values between 1-4
        auto size = []() -> size_t {
            auto ptr = std::getenv("STRINGSIZE");
            if (ptr) {
                return std::stoull(ptr);
            }
            #ifdef NDEBUG
                return 1'000'000;
            #else
                return 1'000;
            #endif
        }();
        return generateText<min, range>(size);
    }();
    return text;
}

template <size_t min, size_t range>
auto generateTexts() -> std::vector<std::vector<uint8_t>> const& {
    static auto text = std::vector<std::vector<uint8_t>>{generateText<min, range>()};
    return text;
}


template <size_t min, size_t range>
auto generateLargeText(size_t length) -> std::vector<uint64_t> {
    auto rng = ankerl::nanobench::Rng{};

    auto text = std::vector<uint64_t>{};
    for (size_t i{0}; i<length; ++i) {
        text.push_back(rng.bounded(range) + min);
    }
    return text;
}

template <size_t min, size_t range>
auto generateLargeText() -> std::vector<uint64_t> const& {
    static auto text = []() -> std::vector<uint64_t> {
        auto rng = ankerl::nanobench::Rng{};

        // generates string with values between 1-4
        auto size = []() -> size_t {
            auto ptr = std::getenv("STRINGSIZE");
            if (ptr) {
                return std::stoull(ptr);
            }
            #ifdef NDEBUG
                return 1'000'000;
            #else
                return 1'000;
            #endif
        }();
        return generateLargeText<min, range>(size);
    }();
    return text;
}

}

#ifdef FMC_USE_AWFMINDEX
    #include "AWFMIndex.h"
    #define STRINGSWITHRANK ALLSTRINGSWITHRANK, AWFMIndex
#else
    #define STRINGSWITHRANK ALLSTRINGSWITHRANK
#endif
