#include "zeroOnesZero.h"

#include <cassert>

namespace search_schemes::generator {

auto zeroOnesZero_trivial(int minK, int K) -> Scheme {
    int N = K+2;

    assert(N>0);
    assert(minK >= 0);
    assert(minK <= K);
    assert(K>=0);
    assert(N>K);

    auto res = Scheme{};

    for (int i{0}; i < N; ++i) {
        for (int j{0}; j < N-i-1; ++j) {
            auto s = Search{};
            // set order
            s.pi.push_back(i);
            for (int k{i+1}; k < N; ++k) {
                s.pi.push_back(k);
            }
            for (int k{i}; k > 0; --k) {
                s.pi.push_back(k-1);
            }

            // set lower bound
            s.l.push_back(0);
            for (int k{0}; k < j; ++k) {
                s.l.push_back(1+k);
            }
            s.l.push_back(j);
            for (int k{j+2}; k < N; ++k) {
                s.l.push_back(j);
            }

            // set upper bound
            s.u.push_back(0);
            for (int k{0}; k < j; ++k) {
                s.u.push_back(1+k);
            }
            s.u.push_back(j);
            for (int k{j+2}; k < N; ++k) {
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

auto zeroOnesZero_opt(int minK, int K) -> Scheme {
    int N = K+2;

    assert(N>0);
    assert(minK >= 0);
    assert(minK <= K);
    assert(K>=0);
    assert(N>K);

    auto res = Scheme{};

    for (int i{0}; i < N-1; ++i) {
        for (int j{0}; j < N-i-1; ++j) {
            auto s = Search{};
            // set order
            s.pi.push_back(i);
            for (int k{i+1}; k < N; ++k) {
                s.pi.push_back(k);
            }
            for (int k{i}; k > 0; --k) {
                s.pi.push_back(k-1);
            }

            // set lower bound
            s.l.push_back(0);
            for (int k{0}; k < j; ++k) {
                s.l.push_back(1+k);
            }
            s.l.push_back(j);
            for (int k{j+2}; k < N; ++k) {
                s.l.push_back(j);
            }

            // set upper bound
            s.u.push_back(0);
            for (int k{0}; k < j; ++k) {
                s.u.push_back(1+k);
            }
            s.u.push_back(j);
            for (int k{j+2}; k < N; ++k) {
                s.u.push_back(K);
            }

            if (not res.empty() and res.back().pi == s.pi) {
                auto& r = res.back();
                for (int k{0}; k < N; ++k) {
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
