// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "FMIndex.h"
#include "BiFMIndex.h"
#include "../string/utils.h"

#include <cassert>
#include <type_traits>
#include <vector>

namespace fmindex_collection {

/**
 * creates the R array for interleaving FM-Indices
 */
template <typename StringLhs, typename StringRhs, typename value_t = size_t>
auto computeInterleavingR(StringLhs const& lhsStr, StringRhs const& rhsStr) -> std::vector<bool> {
    if (rhsStr.size() > std::numeric_limits<value_t>::max()) {
        throw std::runtime_error{"Can not create interleaving R for this value type, index to large"};
    }

    auto lhsC = string::computeAccumulatedC(lhsStr);
    auto rhsC = string::computeAccumulatedC(rhsStr);

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

template <String_c StringLhs, String_c StringRhs>
auto computeMergedBwt(std::vector<bool> const& R, StringLhs const& lhsBwt, StringRhs const& rhsBwt) -> std::vector<uint8_t> {
    auto mergedBwt = std::vector<uint8_t>{};
    mergedBwt.reserve(lhsBwt.size() + rhsBwt.size());

    size_t idx1{}, idx2{};
    for (bool v : R) {
        if (!v) {
            assert(idx1 < lhsBwt.size());
            mergedBwt.push_back(lhsBwt.symbol(idx1));
            idx1 += 1;
        } else {
            assert(idx2 < index2.size());
            mergedBwt.push_back(rhsBwt.symbol(idx2));
            idx2 += 1;
        }
    }
    return mergedBwt;
}

template <typename TCSA>
auto computeCsa(std::vector<bool> const& R, TCSA const& lhsCsa, TCSA const& rhsCsa) -> TCSA {
    auto csa = TCSA::createJoinedCSA(lhsCsa, rhsCsa);
    size_t idx1{}, idx2{};
    for (bool v : R) {
        if (!v) {
            csa.push_back(lhsCsa.value(idx1));
            idx1 += 1;
        } else {
            csa.push_back(rhsCsa.value(idx2));
            idx2 += 1;
        }
    }
    return csa;
}

template <typename Res = void, String_c StringLhs, String_c StringRhs, typename TCSA>
auto mergeImpl(FMIndex<StringLhs, TCSA> const& index1, FMIndex<StringRhs, TCSA> const& index2) -> FMIndex<std::conditional_t<std::is_void_v<Res>, StringLhs, Res>, TCSA> {
    auto R         = computeInterleavingR(index1.bwt, index2.bwt);
    auto mergedBwt = computeMergedBwt(R, index1.bwt, index2.bwt);
    auto csa       = computeCsa(R, index1.csa, index2.csa);

    return {mergedBwt, std::move(csa)};
}

template <typename Res = void, String_c StringLhs, String_c StringRhs, typename TCSA>
auto merge(FMIndex<StringLhs, TCSA> const& index1, FMIndex<StringRhs, TCSA> const& index2) -> FMIndex<std::conditional_t<std::is_void_v<Res>, StringLhs, Res>, TCSA> {
    return mergeImpl<Res>(index1, index2);
}


template <typename Res = void, String_c StrLhs, String_c StrRhs, typename TCSA>
auto mergeImpl(BiFMIndex<StrLhs, TCSA> const& index1, BiFMIndex<StrRhs, TCSA> const& index2) -> BiFMIndex<std::conditional_t<std::is_void_v<Res>, StrLhs, Res>, TCSA> {

    auto R            = computeInterleavingR(index1.bwt, index2.bwt);
    auto mergedBwt    = computeMergedBwt(R, index1.bwt, index2.bwt);
    auto csa          = computeCsa(R, index1.csa, index2.csa);

    R                 = computeInterleavingR(index1.bwtRev, index2.bwtRev);
    auto mergedBwtRev = computeMergedBwt(R, index1.bwtRev, index2.bwtRev);

    return {mergedBwt, mergedBwtRev, std::move(csa)};
}

template <typename Res = void, String_c StrLhs, String_c StrRhs, typename TCSA>
auto merge(BiFMIndex<StrLhs, TCSA> const& index1, BiFMIndex<StrRhs, TCSA> const& index2) -> BiFMIndex<std::conditional_t<std::is_void_v<Res>, StrLhs, Res>, TCSA> {
    return mergeImpl(index1, index2);
}

}
