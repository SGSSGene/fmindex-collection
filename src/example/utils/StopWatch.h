// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <chrono>

struct StopWatch {
    using TP = decltype(std::chrono::steady_clock::now());

    TP start{std::chrono::steady_clock::now()};

    auto now() const -> TP {
        return std::chrono::steady_clock::now();
    }

    auto reset() -> double {
        auto t = now();
        auto diff = t - start;
        start = t;
        return diff.count() / 1'000'000'000.;
    }

    auto peek() const -> double {
        auto t = now();
        auto diff = t - start;
        return diff.count() / 1'000'000'000.;
    }
};

