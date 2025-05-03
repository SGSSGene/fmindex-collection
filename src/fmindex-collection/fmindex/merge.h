// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "FMIndex.h"
#include "BiFMIndex.h"

#include <cassert>
#include <type_traits>
#include <vector>

namespace fmindex_collection {

/**
 * creates the R array for interleaving FM-Indices
 */
template <typename OccLhs, typename OccRhs, typename value_t = size_t>
auto computeInterleavingR(OccLhs const& lhsOcc, OccRhs const& rhsOcc) -> std::vector<bool> {
    if (rhsOcc.size() > std::numeric_limits<value_t>::max()) {
        throw std::runtime_error{"Can not create interleaving R for this value type, index to large"};
    }

    auto R = std::vector<bool>{};
    R.resize(lhsOcc.size() + rhsOcc.size(), false);

    auto nbrOfSeqRhs = rhsOcc.rank(rhsOcc.size(), 0);
    for (size_t n{}; n < nbrOfSeqRhs; ++n) {
        size_t idx1{};
        size_t idx2{n};
        uint8_t c{};
        do {
            assert(idx1 + idx2 < R.size());
            assert(R[idx1 + idx2] == false);
            R[idx1 + idx2] = true;
            c = rhsOcc.symbol(idx2);
            idx1 = lhsOcc.rank(idx1, c);
            idx2 = rhsOcc.rank(idx2, c);
        } while(c != 0);
    }
    assert([&]() {
        size_t a{};
        for (auto b : R) {
            a += b;
        }
        return a;
    }() == rhsOcc.size());
    return R;
}

/**
 * creates the R array for interleaving FM-Indices
 */
template <typename StringLhs, typename StringRhs, typename value_t = size_t>
auto computeInterleavingR(StringLhs const& lhsStr, std::array<size_t, StringLhs::Sigma+1> const& lhsC, StringRhs const& rhsStr, std::array<size_t, StringRhs::Sigma+1> const& rhsC) -> std::vector<bool> {
    if (rhsStr.size() > std::numeric_limits<value_t>::max()) {
        throw std::runtime_error{"Can not create interleaving R for this value type, index to large"};
    }

    auto R = std::vector<bool>{};
    R.resize(lhsStr.size() + rhsStr.size(), false);

    auto nbrOfSeqRhs = rhsStr.rank(rhsStr.size(), 0);
    for (size_t n{}; n < nbrOfSeqRhs; ++n) {
        size_t idx1{};
        size_t idx2{n};
        uint8_t c{};
        do {
            assert(idx1 + idx2 < R.size());
            assert(R[idx1 + idx2] == false);
            R[idx1 + idx2] = true;
            c = rhsStr.symbol(idx2);
            idx1 = lhsStr.rank(idx1, c) + lhsC[c];
            idx2 = rhsStr.rank(idx2, c) + rhsC[c];
        } while(c != 0);
    }
    assert([&]() {
        size_t a{};
        for (auto b : R) {
            a += b;
        }
        return a;
    }() == rhsStr.size());
    return R;
}



template <typename Res = void, typename OccLhs, typename OccRhs, typename TCSA>
auto mergeImpl(FMIndex<OccLhs, TCSA> const& index1, FMIndex<OccRhs, TCSA> const& index2, size_t seqOffset1, size_t seqOffset2) -> FMIndex<std::conditional_t<std::is_void_v<Res>, OccLhs, Res>, TCSA> {
    auto R = computeInterleavingR(index1.bwt, index1.C, index2.bwt, index2.C);

    // Interleave BWT->R and SA->ssa
    auto mergedBWT = std::vector<uint8_t>{};
    mergedBWT.reserve(index1.size() + index2.size());

    using CSA = decltype(index1.csa);
    auto csa = CSA::createJoinedCSA(index1.csa, index2.csa);

    auto addSSAEntry = [&csa](auto const& index, size_t idx, size_t seqOffset) {
        auto loc = index.csa.value(idx);
        if (loc) {
            auto [seq, pos] = *loc;
            csa.push_back(std::tuple{seq + seqOffset, pos});
        } else {
            csa.push_back(std::nullopt);
        }


    };

    size_t idx1{}, idx2{};
    for (bool v : R) {
        if (!v) {
            assert(idx1 < index1.size());
            mergedBWT.push_back(index1.bwt.symbol(idx1));
            addSSAEntry(index1, idx1, seqOffset1);
            idx1 += 1;
        } else {
            assert(idx2 < index2.size());
            mergedBWT.push_back(index2.bwt.symbol(idx2));
            addSSAEntry(index2, idx2, seqOffset2);
            idx2 += 1;
        }
    }
    R.clear();

    return {mergedBWT, std::move(csa)};
}

template <typename Res = void, typename OccLhs, typename OccRhs, typename TCSA>
auto merge(FMIndex<OccLhs, TCSA> const& index1, FMIndex<OccRhs, TCSA> const& index2) -> FMIndex<std::conditional_t<std::is_void_v<Res>, OccLhs, Res>, TCSA> {
//    if (index1.size() >= index2.size()) {
        auto r = index1.bwt.rank(index1.size(), 0) + index1.C[0];
        return mergeImpl<Res>(index1, index2, 0, r);
//    }
//    return mergeImpl<Res>(index2, index1, index2.occ.rank(index2.size(), 0), 0);
}


template <typename Res = void, typename StrLhs, typename StrRhs, typename TCSA>
auto mergeImpl(BiFMIndex<StrLhs, TCSA> const& index1, BiFMIndex<StrRhs, TCSA> const& index2, size_t seqOffset1, size_t seqOffset2) -> BiFMIndex<std::conditional_t<std::is_void_v<Res>, StrLhs, Res>, TCSA> {
    auto csa = TCSA::createJoinedCSA(index1.csa, index2.csa);

    // Interleave BWT->R and SA->ssa
    auto mergedBWT = std::vector<uint8_t>{};
    mergedBWT.reserve(index1.size() + index2.size());

    // compute normal forward bwt
    {
        auto R = computeInterleavingR(index1.bwt, index1.C, index2.bwt, index2.C);
        auto addSSAEntry = [&csa](auto const& index, size_t idx, size_t seqOffset) {
            auto loc = index.csa.value(idx);
            if (loc) {
                auto [seq, pos] = *loc;
                csa.push_back(std::tuple{seq + seqOffset, pos});
            } else {
                csa.push_back(std::nullopt);
            }
        };
        size_t idx1{}, idx2{};
        for (bool v : R) {
            if (!v) {
                assert(idx1 < index1.bwt.size());
                mergedBWT.push_back(index1.bwt.symbol(idx1));
                addSSAEntry(index1, idx1, seqOffset1);
                idx1 += 1;
            } else {
                assert(idx2 < index2.bwt.size());
                mergedBWT.push_back(index2.bwt.symbol(idx2));
                addSSAEntry(index2, idx2, seqOffset2);
                idx2 += 1;
            }
        }
    }

    // Interleave BWT->R and SA->ssa
    auto mergedBWTRev = std::vector<uint8_t>{};
    mergedBWTRev.reserve(mergedBWT.size());

    // compute reversed bwt
    {
        auto R = computeInterleavingR(index1.bwtRev, index1.C, index2.bwtRev, index2.C);
        size_t idx1{}, idx2{};
        for (bool v : R) {
            if (!v) {
                assert(idx1 < index1.bwtRev.size());
                mergedBWTRev.push_back(index1.bwtRev.symbol(idx1));
                idx1 += 1;
            } else {
                assert(idx2 < index2.bwtRev.size());
                mergedBWTRev.push_back(index2.bwtRev.symbol(idx2));
                idx2 += 1;
            }
        }
    }

    return {mergedBWT, mergedBWTRev, std::move(csa)};
}

template <typename Res = void, typename StrLhs, typename StrRhs, typename TCSA>
auto merge(BiFMIndex<StrLhs, TCSA> const& index1, BiFMIndex<StrRhs, TCSA> const& index2) -> BiFMIndex<std::conditional_t<std::is_void_v<Res>, StrLhs, Res>, TCSA> {
//    if (index1.size() >= index2.size()) {
        return mergeImpl(index1, index2, 0, index1.bwt.rank(index1.size(), 0) + index1.C[0]);
//    }
//    return mergeImpl(index2, index1, index2.occ.rank(index2.size(), 0), 0);
}

}
