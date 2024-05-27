// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "FMIndex.h"
#include "BiFMIndex.h"

#include <type_traits>
#include <vector>

namespace fmindex_collection {

/**
 * creates the R array for interleaving FMIndices
 */
template <typename OccLhs, typename OccRhs, typename TCSA, typename value_t = size_t>
auto computeInterleavingR(FMIndex<OccLhs, TCSA> const& index1, FMIndex<OccRhs, TCSA> const& index2) -> std::vector<value_t> {
    auto R = std::vector<value_t>{};
    R.resize(index2.size(), 0);

    if (index2.size() > std::numeric_limits<value_t>::max()) {
        throw std::runtime_error{"Can not create interleaving R for this value type, index to large"};
    }

    size_t idx1{};
    size_t idx2{};

    for (size_t i{0}; i < R.size(); ++i) {
        auto c = index2.occ.symbol(idx2);
        idx1 = index1.occ.rank(idx1, c);
        idx2 = index2.occ.rank(idx2, c);
        R[idx2] = idx1;
    }
    return R;
}



//!TODO swap index1 and index2 if index1 is smaller
template <typename Res = void, typename OccLhs, typename OccRhs, typename TCSA>
auto merge(FMIndex<OccLhs, TCSA> const& index1, FMIndex<OccRhs, TCSA> const& index2) -> FMIndex<std::conditional_t<std::is_void_v<Res>, OccLhs, Res>, TCSA> {

    auto R = computeInterleavingR(index1, index2);

    // Interleave BWT->R and SA->ssa
    auto mergedBWT = std::vector<uint8_t>{};
    mergedBWT.reserve(index1.size() + index2.size());

    // The right hand index sequences, require adjusted sequence number
    size_t seqOffset = index1.occ.rank(index1.size(), 0);

    using CSA = decltype(index1.csa);
    auto csa = CSA::createJoinedCSA(index1.csa, index2.csa);

    auto addSSAEntry = [&csa](auto const& index, size_t idx, size_t seqOffset) {
        auto loc = index.csa.value(idx);
        if (loc) {
            auto [seq, pos] = *loc;
            csa.push_back(std::tuple{seq + seqOffset, pos});
        }
    };

    size_t idx1{};
    for (size_t idx2{}; idx2 < R.size(); ++idx2) {
        for (; idx1 < R[idx2]; ++idx1) {
            mergedBWT.push_back(index1.occ.symbol(idx1));
            addSSAEntry(index1, idx1, 0);
        }
        mergedBWT.push_back(index2.occ.symbol(idx2));
        addSSAEntry(index2, idx2, seqOffset);
    }
    for (; idx1 < index1.size(); ++idx1) {
        mergedBWT.push_back(index1.occ.symbol(idx1));
        addSSAEntry(index1, idx1, 0);
    }


    return {mergedBWT, std::move(csa)};
}

}
