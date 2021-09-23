#pragma once

#include <tuple>
#include <vector>

namespace search_schemes {

struct SearchTree {
    std::vector<int> pi;
    std::vector<int> l;
    std::vector<int> u;

    bool operator==(SearchTree const& _other) const {
        return std::tie(pi, l, u) == std::tie(_other.pi, _other.l, _other.u);
    }
};

}
