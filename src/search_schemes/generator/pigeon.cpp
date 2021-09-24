#include "pigeon.h"

#include <cassert>

namespace search_schemes::generator {

auto pigeon_trivial(int minK, int K) -> Scheme {
    int N = K+1;

    assert(N>0);
    assert(minK >= 0);
    assert(minK <= K);
    assert(K>=0);
    assert(N>K);

    auto res = Scheme{};

    for (int i{0}; i < N; ++i) {
        auto s = Search{};
        // set order
        s.pi.push_back(i);
        for (int j{i-1}; j >= 0; --j) {
            s.pi.push_back(j);
        }
        for (int j{i+1}; j < N; ++j) {
            s.pi.push_back(j);
        }


        // set lower bound
        for (int j{0}; j < N; ++j) {
            s.l.push_back(0);
        }

        // set upper bound
        s.u.push_back(0);
        for (int j{1}; j < N; ++j) {
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

auto pigeon_opt(int minK, int K) -> Scheme {
    int N = K+1;

    assert(N>0);
    assert(minK >= 0);
    assert(minK <= K);
    assert(K>=0);
    assert(N>K);

    auto res = Scheme{};

    for (int i{0}; i < N; ++i) {
        auto s = Search{};
        // set order
        s.pi.push_back(i);
        for (int j{i-1}; j >= 0; --j) {
            s.pi.push_back(j);
        }
        for (int j{i+1}; j < N; ++j) {
            s.pi.push_back(j);
        }


        // set lower bound
        s.l.push_back(0);
        for (int j{i-1}; j >= 0; --j) {
            s.l.push_back(i-j);
        }
        for (int j{i+1}; j < N; ++j) {
            s.l.push_back(i);
        }

        // set upper bound
        s.u.push_back(0);
        for (int j{i-1}; j >= 0; --j) {
            s.u.push_back(std::min(K, K-j));
        }
        for (int j{i+1}; j < N; ++j) {
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
