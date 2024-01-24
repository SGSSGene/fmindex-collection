// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../concepts.h"
#include "FMIndex.h"
#include "FMIndexCursor.h"
#include "occtable/Bitvector.h"

namespace fmindex_collection::search_no_errors_fpga {

using index_t = fmindex_collection::FMIndex<occtable::bitvector::OccTable<5>>;

/* query only has values between 1-4
 *
 * returns the range of valid values
 */
inline auto search(index_t const & index, std::string const& query) -> std::tuple<size_t, size_t> {

    size_t lb = 0;
    size_t rb = index.size(); // rb is exclusive boundary

    for (size_t i{0}; i < query.size(); ++i) {
        auto symb = query[query.size() - i - 1];
        { // extend to the left
            size_t newLb  = index.occ.rank(lb, symb);
            size_t newLen = index.occ.rank(rb, symb) - newLb;
            size_t newRb = newLb + newLen;
            lb = newLb;
            rb = newRb;
        }
        if (lb == rb) { // check if its empty
            break;
        }
    }
    return {lb, rb};
}

}
