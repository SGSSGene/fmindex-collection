// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#pragma once

#include <catch2/catch_all.hpp>
#include <fmindex-collection/utils.h>
#include <fstream>
#include <nanobench.h>

#include "../BenchSize.h"
#include "allStrings.h"

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

template <typename T>
struct X;
template <typename ...T>
struct X<std::tuple<T...>> {
    void operator()(auto f) {
        (f.template operator()<T>(), ...);
    }
};

template <typename Tuple>
void call_with_templates_as_tuple(auto f) {
    X<Tuple>{}(f);
}

template <typename ...T>
void call_with_templates(auto f) {
    call_with_templates_as_tuple<std::tuple<T...>>(f);
}

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
            auto ptr = std::getenv("VECTORSIZE");
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
}
