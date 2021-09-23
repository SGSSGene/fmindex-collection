#include "isValid.h"

#include <cassert>
#include <algorithm>

namespace oss {

namespace {
    // check if pi is contiguous and start with 1
    auto checkPiContiguous(decltype(SearchTree::pi) pi) {
        assert(not pi.empty());
        std::sort(begin(pi), end(pi));
        auto iter = std::unique(begin(pi), end(pi));
        if (iter != end(pi)) {
            return false;
        }
        if (pi.front() != 1) {
            return false;
        }
        if (static_cast<size_t>(pi.back()) != pi.size()) {
            return false;
        }
        return true;
    }

    // check that if is monotonically increasing
    template <typename Collection>
    auto checkMonotonicallyIncreasing(Collection const& col) {
        assert(not col.empty());
        if (col.size() == 1) {
            return true;
        }
        auto iter = begin(col);
        for (auto lastIter = iter++; iter != end(col); lastIter = iter++) {
            if (*lastIter > *iter) {
                return false;
            }
        }
        return true;
    }

}

auto isValid(SearchTree const& s) -> bool {
    // check if has valid size
    if (s.pi.empty()) {
        return false;
    }

    // check if all entries have same length
    if (s.pi.size() != s.u.size() or s.pi.size() != s.l.size()) {
        return false;
    }

    // check if pi is contiguous and start with 1
    if (not checkPiContiguous(s.pi)) {
        return false;
    }
    // check that l and u are monotonically increasing
    if (not checkMonotonicallyIncreasing(s.l)) {
        return false;
    }
    if (not checkMonotonicallyIncreasing(s.u)) {
        return false;
    }
    // check that l is always equal or smaller than u
    for (size_t i{0}; i < s.pi.size(); ++i) {
        if (s.l[i] > s.u[i]) {
            return false;
        }
    }
    return true;
}

auto isValid(Scheme const& ss) -> bool {
    // searching for an invalid scheme
    auto iter = std::find_if(begin(ss), end(ss), [&](auto const& s) {
        return not isValid(s) or s.pi.size() != ss.front().pi.size();
    });
    return iter == end(ss);
}

}
