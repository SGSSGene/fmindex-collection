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
template <typename OccLhs, typename OccRhs, typename value_t = size_t>
auto computeInterleavingR(OccLhs const& lhsOcc, OccRhs const& rhsOcc) -> std::vector<value_t> {
    if (rhsOcc.size() > std::numeric_limits<value_t>::max()) {
        throw std::runtime_error{"Can not create interleaving R for this value type, index to large"};
    }

    auto R = std::vector<value_t>{};
    R.resize(rhsOcc.size(), 0);

    size_t idx1{};
    size_t idx2{};

    for (size_t i{0}; i < R.size(); ++i) {
        auto c = rhsOcc.symbol(idx2);
        idx1 = lhsOcc.rank(idx1, c);
        idx2 = rhsOcc.rank(idx2, c);
        R[idx2] = idx1;
    }
    return R;
}



//!TODO swap index1 and index2 if index1 is smaller
template <typename Res = void, typename OccLhs, typename OccRhs, typename TCSA>
auto merge(FMIndex<OccLhs, TCSA> const& index1, FMIndex<OccRhs, TCSA> const& index2) -> FMIndex<std::conditional_t<std::is_void_v<Res>, OccLhs, Res>, TCSA> {
    auto R = computeInterleavingR(index1.occ, index2.occ);

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
    R.clear();


    return {mergedBWT, std::move(csa)};
}

//!TODO swap index1 and index2 if index1 is smaller
template <typename Res = void, typename OccLhs, typename OccRhs, typename TCSA>
auto mergeImpl(BiFMIndex<OccLhs, TCSA> const& index1, BiFMIndex<OccRhs, TCSA> const& index2, size_t seqOffset1, size_t seqOffset2) -> BiFMIndex<std::conditional_t<std::is_void_v<Res>, OccLhs, Res>, TCSA> {
    auto csa = TCSA::createJoinedCSA(index1.csa, index2.csa);

    // Interleave BWT->R and SA->ssa
    auto mergedBWT = std::vector<uint8_t>{};
    mergedBWT.reserve(index1.size() + index2.size());

    // compute normal forward bwt
    {
        auto R = computeInterleavingR(index1.occ, index2.occ);
/*
        // The right hand index sequences, require adjusted sequence number
        auto [seqOffset1, seqOffset2] = [&]() -> std::tuple<size_t, size_t> {
            if (swapOffsets) {
                return {index2.occ.rank(index2.size(), 0), 0};
            }
            return {0, index1.occ.rank(index1.size(), 0)};
        }();
*/

        auto addSSAEntry = [&csa](auto const& index, size_t idx, size_t seqOffset) {
            auto loc = index.csa.value(idx);
            if (loc) {
                auto [seq, pos] = *loc;
                csa.push_back(std::tuple{seq + seqOffset, pos});
            } else {
                csa.push_back(std::nullopt);
            }
        };

        size_t idx1{};
        for (size_t idx2{}; idx2 < R.size(); ++idx2) {
            for (; idx1 < R[idx2]; ++idx1) {
                mergedBWT.push_back(index1.occ.symbol(idx1));
                addSSAEntry(index1, idx1, seqOffset1);
            }
            mergedBWT.push_back(index2.occ.symbol(idx2));
            addSSAEntry(index2, idx2, seqOffset2);
        }
        for (; idx1 < index1.size(); ++idx1) {
            mergedBWT.push_back(index1.occ.symbol(idx1));
            addSSAEntry(index1, idx1, seqOffset1);
        }
    }

    // Interleave BWT->R and SA->ssa
    auto mergedBWTRev = std::vector<uint8_t>{};
    mergedBWTRev.reserve(index1.size() + index2.size());

    // compute reversed bwt
    {
        auto R = computeInterleavingR(index1.occRev, index2.occRev);

        size_t idx1{};
        for (size_t idx2{}; idx2 < R.size(); ++idx2) {
            for (; idx1 < R[idx2]; ++idx1) {
                mergedBWTRev.push_back(index1.occRev.symbol(idx1));
            }
            mergedBWTRev.push_back(index2.occRev.symbol(idx2));
        }
        for (; idx1 < index1.size(); ++idx1) {
            mergedBWTRev.push_back(index1.occRev.symbol(idx1));
        }
    }

    return {mergedBWT, mergedBWTRev, std::move(csa)};
}

template <typename Res = void, typename OccLhs, typename OccRhs, typename TCSA>
auto merge(BiFMIndex<OccLhs, TCSA> const& index1, BiFMIndex<OccRhs, TCSA> const& index2) -> BiFMIndex<std::conditional_t<std::is_void_v<Res>, OccLhs, Res>, TCSA> {
    if (index1.size() >= index2.size()) {
        return mergeImpl(index1, index2, 0, index1.occ.rank(index1.size(), 0));
    }
    return mergeImpl(index2, index1, index2.occ.rank(index2.size(), 0), 0);
}

}
