#pragma once

#include <tuple>
#include <vector>

namespace search_schemes {

/* Single Search
 *
 * Each search is defined by pi indicating the order of the parts.
 * Each part has an individual lower and upper limit l and u.
 *
 * Note: pi.size() == l.size() == u.size()
 */
struct Search {
    std::vector<int> pi;
    std::vector<int> l;
    std::vector<int> u;

    bool operator==(Search const& _other) const {
        return std::tie(pi, l, u) == std::tie(_other.pi, _other.l, _other.u);
    }
};

}
