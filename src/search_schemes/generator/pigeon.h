// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../Scheme.h"

#include <cassert>

namespace search_schemes::generator {

inline auto pigeon_trivial(size_t minK, size_t K) -> Scheme {
    auto N = K+1;

    assert(N>0);
    assert(minK <= K);
    assert(N>K);

    auto res = Scheme{};

    for (size_t i{0}; i < N; ++i) {
        auto s = Search{};
        // set order
        s.pi.push_back(i);
        for (size_t j{i}; j > 0; --j) {
            s.pi.push_back(j-1);
        }
        for (size_t j{i+1}; j < N; ++j) {
            s.pi.push_back(j);
        }


        // set lower bound
        for (size_t j{0}; j < N; ++j) {
            s.l.push_back(0);
        }

        // set upper bound
        s.u.push_back(0);
        for (size_t j{1}; j < N; ++j) {
            s.u.push_back(K);
        }
        res.push_back(s);
    }

    // set minimum value
    for (auto& s : res) {
        s.l.back() = std::max(s.l.back(), minK);
    }

    return res;
}

inline auto pigeon_opt(size_t minK, size_t K) -> Scheme {
    auto N = K+1;

    assert(N>0);
    assert(minK <= K);
    assert(N>K);

    auto res = Scheme{};

    for (size_t i{0}; i < N; ++i) {
        auto s = Search{};
        // set order
        s.pi.push_back(i);
        for (size_t j{i}; j > 0; --j) {
            s.pi.push_back(j-1);
        }
        for (size_t j{i+1}; j < N; ++j) {
            s.pi.push_back(j);
        }


        // set lower bound
        s.l.push_back(0);
        for (size_t j{i}; j > 0; --j) {
            s.l.push_back(i-j+1);
        }
        for (size_t j{i+1}; j < N; ++j) {
            s.l.push_back(i);
        }

        // set upper bound
        s.u.push_back(0);
        for (size_t j{i}; j > 0; --j) {
            s.u.push_back(K-j+1);
        }
        for (size_t j{i+1}; j < N; ++j) {
            s.u.push_back(K);
        }

        res.push_back(s);
    }

    // set minimum value
    for (auto& s : res) {
        s.l.back() = std::max(s.l.back(), minK);
    }

    return res;
}
}
