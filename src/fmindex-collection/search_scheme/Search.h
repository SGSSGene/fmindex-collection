// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <cstddef>
#include <tuple>
#include <vector>

namespace fmc::search_scheme {

/* Single Search
 *
 * Each search is defined by pi indicating the order of the parts.
 * Each part has an individual lower and upper limit l and u.
 *
 * Note: pi.size() == l.size() == u.size()
 */
struct Search {
    std::vector<size_t> pi;
    std::vector<size_t> l;
    std::vector<size_t> u;

    bool operator==(Search const& _other) const {
        return std::tie(pi, l, u) == std::tie(_other.pi, _other.l, _other.u);
    }
};

}
