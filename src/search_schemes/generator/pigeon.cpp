#include "pigeon.h"

#include "../nodeCount.h"

#include <cassert>
#include <limits>

namespace oss::generator {

namespace {

auto orderAndAccumulateByProjection(SearchTree s, int K) {
    auto l = s.l;
    auto u = s.u;
    s.l.clear();
    s.u.clear();
    auto inverse_pi = [&](int i) {
        for (int j{0}; j<s.pi.size(); ++j) {
            if (s.pi[j] == i) {
                return j;
            }
        }
        assert(false);
    };

    // reorder upper and lower bound and accumulate
    int acc_l{}, acc_u{};
    for (int p{0}; p < s.pi.size(); ++p) {
        acc_l = std::min(acc_l + l[s.pi[p]], K);
        acc_u = std::min(acc_u + u[s.pi[p]], K);
        s.l.push_back(acc_l);
        s.u.push_back(acc_u);
    }
    return s;
}


auto findCost(SearchTree s, int K, int startPos, int sigma) {
    s.pi.clear();
    auto s1{s};
    auto s2{s};

    // variant 1
    for (int i{startPos}; i < s.l.size(); ++i) {
        s1.pi.push_back(i);
    }
    for (int i{startPos-1}; i >= 0; --i) {
        s1.pi.push_back(i);
    }

    // variant 2
    for (int i{startPos}; i >= 0; --i) {
        s2.pi.push_back(i);
    }
    for (int i{startPos+1}; i < s.l.size(); ++i) {
        s2.pi.push_back(i);
    }


    s1 = orderAndAccumulateByProjection(s1, K);
    s2 = orderAndAccumulateByProjection(s2, K);

    for (auto& v : s1.pi) {
        v += 1;
    }
    for (auto& v : s2.pi) {
        v += 1;
    }


    int c1 = nodeCount(std::vector{s1}, sigma);
    if (c1 == 0) {
        c1 = std::numeric_limits<int>::max();
    }
    int c2 = nodeCount(std::vector{s2}, sigma);
    if (c2 == 0) {
        c2 = std::numeric_limits<int>::max();
    }
    if (c2 < c1) {
        return std::make_tuple(c2, s2);
    }
    return std::make_tuple(c1, s1);
}


auto findMinimalCost(SearchTree const& s, int K, int sigma) {
    auto res = SearchTree{};
    int min_c = std::numeric_limits<int>::max();
    for (int i{0}; i < s.l.size(); ++i) {
        auto [c, order] = findCost(s, K, i, sigma);
        if (c < min_c) {
            min_c = c;
            res = order;
        }
    }
    return res;
}
}

auto pigeon(int minK, int K, int sigma) -> Scheme {
    int N = K+1;

    assert(N>0);
    assert(minK >= 0);
    assert(minK <= K);
    assert(K>=0);
    assert(N>K);

    auto res = Scheme{};

    for (int i{0}; i < N; ++i) {
        auto s = SearchTree{};

        // set lower bound
        for (int j{0}; j < i; ++j) {
            s.l.push_back(1);
        }
        for (int j{i}; j < N; ++j) {
            s.l.push_back(0);
        }
        // set upper bound
        for (int j{0}; j < i; ++j) {
            s.u.push_back(K-i+1);
        }
        s.u.push_back(0);
        for (int j{i+1}; j < N; ++j) {
            s.u.push_back(K-i);
        }
        s = findMinimalCost(s, K, sigma);
        res.push_back(s);
    }

    // set minimum value
    for (auto& s : res) {
        s.l.back() = std::max(s.l.back(), minK);
    }

    return res;
}
auto pigeon_trivial(int minK, int K) -> Scheme {
    int N = K+1;

    assert(N>0);
    assert(minK >= 0);
    assert(minK <= K);
    assert(K>=0);
    assert(N>K);

    auto res = Scheme{};

    for (int i{0}; i < N; ++i) {
        auto s = SearchTree{};
        // set order
        s.pi.push_back(i+1);
        for (int j{i-1}; j >= 0; --j) {
            s.pi.push_back(j+1);
        }
        for (int j{i+1}; j < N; ++j) {
            s.pi.push_back(j+1);
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
        auto s = SearchTree{};
        // set order
        s.pi.push_back(i+1);
        for (int j{i-1}; j >= 0; --j) {
            s.pi.push_back(j+1);
        }
        for (int j{i+1}; j < N; ++j) {
            s.pi.push_back(j+1);
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
