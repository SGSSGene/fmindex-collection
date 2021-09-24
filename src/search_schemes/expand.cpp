#include "expand.h"

#include "isValid.h"

namespace search_schemes {
namespace {

auto forwards(std::vector<size_t> const& pi) -> std::vector<size_t> {
    auto rets = std::vector<size_t>{};
    rets.push_back(pi.size() == 1 or pi[1] > pi[0]?1:0);
    for (size_t i{1}; i < pi.size(); ++i) {
        rets.push_back(pi[i] > pi[i-1]?1:0);
    }
    return rets;
}

}

auto expandCount(size_t oldLen, size_t newLen) -> std::vector<size_t> {
    assert(oldLen > 0);
    assert(newLen > 0);

    auto const block_size = newLen / oldLen;
    auto const rest       = newLen % oldLen;
    auto counts = std::vector<size_t>(oldLen, block_size);
    for (size_t i{0}; i < rest; ++i) {
        counts[i] += 1;
    }
    return counts;
}

bool isExpandable(Search s, size_t newLen) {
    assert(isValid(s));
    auto counts = expandCount(s.pi.size(), newLen);
    // check that all parts are still satisfiable

    if (s.l[0] > counts[s.pi[0]]) {
        return false;
    }
    for (size_t i{1}; i<s.pi.size(); ++i) {
        if (s.l[i]-s.u[i] > counts[s.pi[i]]) {
            return false;
        }
    }
    return true;
}

auto expandPI(std::vector<size_t> const& pi, size_t _newLen) -> std::vector<size_t> {
    auto counts = expandCount(pi.size(), _newLen);
    auto starts = std::vector<size_t>(pi.size(), 1);
    for (size_t i{1}; i < pi.size(); ++i) {
        starts[i] = starts[i-1] + counts[i-1];
    }
    size_t newLen = std::accumulate(begin(counts), end(counts), 0);

    auto nums = std::vector<size_t>{};
    auto expandForwards = [&](size_t l, size_t u) {
        for (size_t j = l; j <= u; ++j) {
            nums.push_back(j);
        }
    };
    auto expandBackwards = [&](size_t l, size_t u) {
        for (size_t j{u+1}; j > l; --j) {
            nums.push_back(j-1);
        }
    };
    auto expand = [&](size_t i, bool forward) {
        auto l = starts[pi[i]-1];
        auto u = l + counts[pi[i]-1]-1;
        if (forward) expandForwards(l, u);
        else expandBackwards(l, u);
    };

    auto fs = forwards(pi);

    for (size_t i{0}; i < pi.size(); ++i) {
        expand(i, fs[i]);
    }

    return nums;
}


/*auto expandLowerBound(std::vector<int> const& pi, std::vector<int> bound, int _newLen) -> std::vector<int> {
    auto counts = expandCount(pi.size(), _newLen);
    auto starts = std::vector<int>(pi.size(), 1);
    for (size_t i{1}; i < pi.size(); ++i) {
        starts[i] = starts[i-1] + counts[i-1];
    }
    int newLen = std::accumulate(begin(counts), end(counts), 0);
    auto expandedBound = std::vector<int>(newLen, bound.back());

    int j{0};
    for (int i{0}; i < pi.size(); ++i) {
        auto count = counts[pi[i]-1];
        for(;count > 0; --count, ++j) {
            expandedBound[j] = bound[i] - count+1;
        }
    }
    expandedBound.back() = bound.back();

    for (size_t i{1}; i < expandedBound.size(); ++i) {
        expandedBound[i] = std::max(expandedBound[i-1], expandedBound[i]);
    }
    for (size_t i{expandedBound.size()-1}; i > 0; --i) {
        expandedBound[i-1] = std::max(0, std::max(expandedBound[i-1], expandedBound[i]-1));
    }

    return expandedBound;
}

auto expandUpperBound(std::vector<int> const& pi, std::vector<int> bound, int _newLen) -> std::vector<int> {
    auto counts = expandCount(pi.size(), _newLen);
    auto starts = std::vector<int>(pi.size(), 1);
    for (size_t i{1}; i < pi.size(); ++i) {
        starts[i] = starts[i-1] + counts[i-1];
    }
    int newLen = std::accumulate(begin(counts), end(counts), 0);
    auto expandedBound = std::vector<int>(newLen, bound.back());

    int j{0};
    for (int i{0}; i < pi.size(); ++i) {
        auto count = counts[pi[i]-1];
        for(;count > 0; --count, ++j) {
            expandedBound[j] = bound[i];
        }
    }

    expandedBound[0] = std::min(expandedBound[0], 1);
    for (size_t i{1}; i < expandedBound.size(); ++i) {
        expandedBound[i] = std::min(expandedBound[i-1]+1, expandedBound[i]);
    }
    return expandedBound;
}*/

auto expandLowerBound(std::vector<size_t> const& pi, std::vector<size_t> bound, size_t _newLen) -> std::vector<size_t> {
    auto expandedBound = std::vector<size_t>{};

    auto counts = expandCount(pi.size(), _newLen);
    for (size_t i{0}; i < pi.size(); ++i) {
        auto count = counts[pi[i]-1];
        while(count > 1) {
            --count;
            expandedBound.push_back((i>0)?bound[i-1]:0);
        }
        if (count > 0) {
            expandedBound.push_back(bound[i]);
        } else if (count == 0 and not expandedBound.empty()) {
            expandedBound.back() = bound[i];
        }
    }
    return expandedBound;
}

auto expandUpperBound(std::vector<size_t> const& pi, std::vector<size_t> bound, size_t _newLen) -> std::vector<size_t> {
    auto expandedBound = std::vector<size_t>{};

    auto counts = expandCount(pi.size(), _newLen);
    for (size_t i{0}; i < pi.size(); ++i) {
        auto count = counts[pi[i]-1];
        while (count --> 0) {
            expandedBound.push_back(bound[i]);
        }
    }
    return expandedBound;
}


auto expand(Search s, size_t newLen) -> std::optional<Search> {
    auto r = Search{};
    r.pi = expandPI(s.pi, newLen);
    r.l  = expandLowerBound(s.pi, s.l, newLen);
    r.u  = expandUpperBound(s.pi, s.u, newLen);
    if (not isValid(r)) {
        return std::nullopt;
    }
    return {r};
}
auto expand(Scheme ss, size_t newLen) -> Scheme {
    auto r = Scheme{};
    for (auto const& s : ss) {
        auto o = expand(s, newLen);
        if (o) {
            r.push_back(*o);
        }
    }
    return r;
}



}
