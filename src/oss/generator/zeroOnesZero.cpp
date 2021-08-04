#include "zeroOnesZero.h"

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


auto findCost(SearchTree s, int K, int startPos, int sigma, int N) {
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


    auto c1 = nodeCount(std::vector{s1}, sigma);
    assert(c1 > 0);
    if (c1 == 0) {
        c1 = std::numeric_limits<decltype(c1)>::max();
    }
    auto c2 = nodeCount(std::vector{s2}, sigma);
    assert(c2 > 0);

    if (c2 == 0) {
        c2 = std::numeric_limits<decltype(c2)>::max();
    }
    if (c2 < c1) {
        return std::make_tuple(c2, s2);
    }
    return std::make_tuple(c1, s1);
}


auto findMinimalCost(SearchTree const& s, int K, int sigma, int N) {
    auto res = SearchTree{};
    double min_c = std::numeric_limits<double>::max();
    for (int i{0}; i < s.l.size(); ++i) {
        auto [c, order] = findCost(s, K, i, sigma, N);
        if (c < min_c) {
            min_c = c;
            res = order;
        }
    }
    return res;
}
}
auto zeroOnesZero(int minK, int K, int sigma, int DB_SIZE) -> Scheme {
    assert(minK >= 0);
    assert(minK <= K);
    assert(K>=0);

    int N = K+2;
    auto res = Scheme{};

    for (int i{0}; i < N; ++i) {
        for (int j{i+1}; j < N; ++j) {

            std::vector<int> l(N, 0), u(N, 0);
            int includedErrors{};
            for (int k{0}; k < i; ++k) {
                l[k] = 0;
            }
            if (i>0) {
                l[i-1] = 1;
                u[i-1] = 1;
                includedErrors += 1;
            }
            for (int k{i+1}; k < j; ++k) {
                l[k] = 1;
                u[k] = 1;
                includedErrors += 1;
            }
            for (int k{0}; k < i; ++k) {
                u[k] += K - includedErrors;
            }
            for (int k{j+1}; k < N; ++k) {
                u[k] += K - includedErrors;
            }
            auto s = SearchTree{};
            s.l = l;
            s.u = u;
            s = findMinimalCost(s, K, sigma, DB_SIZE);

            res.push_back(s);
        }
    }

    // set minimum value
    for (auto& s : res) {
        s.l.back() = std::max(s.l.back(), minK);
    }

    return res;
}

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
            auto s = SearchTree{};
            // set order
            s.pi.push_back(i+1);
            for (int k{i+1}; k < N; ++k) {
                s.pi.push_back(k+1);
            }
            for (int k{i}; k > 0; --k) {
                s.pi.push_back(k);
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
            auto s = SearchTree{};
            // set order
            s.pi.push_back(i+1);
            for (int k{i+1}; k < N; ++k) {
                s.pi.push_back(k+1);
            }
            for (int k{i}; k > 0; --k) {
                s.pi.push_back(k);
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
