// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "../Scheme.h"

#include <cassert>

namespace search_schemes::generator {

inline auto zeroOnesZero_trivial(size_t minK, size_t K) -> Scheme {
    auto N = K+2;

    assert(N>0);
    assert(minK <= K);
    assert(N>K);

    auto res = Scheme{};

    for (size_t i{0}; i < N; ++i) {
        for (size_t j{0}; j < N-i-1; ++j) {
            auto s = Search{};
            // set order
            s.pi.push_back(i);
            for (size_t k{i+1}; k < N; ++k) {
                s.pi.push_back(k);
            }
            for (size_t k{i}; k > 0; --k) {
                s.pi.push_back(k-1);
            }

            // set lower bound
            s.l.push_back(0);
            for (size_t k{0}; k < j; ++k) {
                s.l.push_back(1+k);
            }
            s.l.push_back(j);
            for (size_t k{j+2}; k < N; ++k) {
                s.l.push_back(j);
            }

            // set upper bound
            s.u.push_back(0);
            for (size_t k{0}; k < j; ++k) {
                s.u.push_back(1+k);
            }
            s.u.push_back(j);
            for (size_t k{j+2}; k < N; ++k) {
                s.u.push_back(K);
            }
            res.push_back(s);
        }
    }

    // set minimum value
    for (auto& s : res) {
        s.l.back() = std::max(s.l.back(), minK);
    }

    return res;
}

inline auto zeroOnesZero_opt(size_t minK, size_t K) -> Scheme {
    auto N = K+2;

    assert(N>0);
    assert(minK <= K);
    assert(N>K);

    auto res = Scheme{};

    for (size_t i{0}; i < N-1; ++i) {
        for (size_t j{0}; j < N-i-1; ++j) {
            auto s = Search{};
            // set order
            s.pi.push_back(i);
            for (size_t k{i+1}; k < N; ++k) {
                s.pi.push_back(k);
            }
            for (size_t k{i}; k > 0; --k) {
                s.pi.push_back(k-1);
            }

            // set lower bound
            s.l.push_back(0);
            for (size_t k{0}; k < j; ++k) {
                s.l.push_back(1+k);
            }
            s.l.push_back(j);
            for (size_t k{j+2}; k < N; ++k) {
                s.l.push_back(j);
            }

            // set upper bound
            s.u.push_back(0);
            for (size_t k{0}; k < j; ++k) {
                s.u.push_back(1+k);
            }
            s.u.push_back(j);
            for (size_t k{j+2}; k < N; ++k) {
                s.u.push_back(K);
            }

            if (not res.empty() and res.back().pi == s.pi) {
                auto& r = res.back();
                for (size_t k{0}; k < N; ++k) {
                    r.l[k] = std::min(r.l[k], s.l[k]);
                    r.u[k] = std::max(r.u[k], s.u[k]);
                }
            } else {
                res.push_back(s);
            }
        }
    }

    // set minimum value
    for (auto& s : res) {
        s.l.back() = std::max(s.l.back(), minK);
    }

    return res;
}

}
