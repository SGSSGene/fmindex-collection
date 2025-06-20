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

namespace fmindex_collection::fmindex {

namespace detail {

/**
 * Creates the R array for interleaving SA/BWT/FM-Indices.
 * \param lhsStr first bwt
 * \param rhsStr second bwt
 * \return an boolean array indicating in which slot which entry is expected
 *
 * The R array has the same size as lhsStr and rhsStr together. It indicates
 * with 'false' that an entry from lhsStr is expected and with 'true' an entry from rhsStr.
 *
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

/** Creates a new BWT from two bwt that are merged with the help of an R array
 */
template <String_c StringLhs, String_c StringRhs>
auto mergeBwt(std::vector<bool> const& R, StringLhs const& lhsBwt, StringRhs const& rhsBwt) -> std::vector<uint8_t> {
    auto mergedBwt = std::vector<uint8_t>{};
    mergedBwt.reserve(lhsBwt.size() + rhsBwt.size());

    size_t idx1{}, idx2{};
    for (bool v : R) {
        if (!v) {
            assert(idx1 < lhsBwt.size());
            mergedBwt.push_back(lhsBwt.symbol(idx1));
            idx1 += 1;
        } else {
            assert(idx2 < rhsBwt.size());
            mergedBwt.push_back(rhsBwt.symbol(idx2));
            idx2 += 1;
        }
    }
    return mergedBwt;
}

/** Creates a new CSA from two csa that are merged with the help of an R array
 */
template <typename TCSA>
auto mergeCsa(std::vector<bool> const& R, TCSA const& lhsCsa, TCSA const& rhsCsa) -> TCSA {
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

template <typename T>
auto mergeSparseArrays(std::vector<bool> const& R, suffixarray::SparseArray<T> const& lhs, suffixarray::SparseArray<T> const& rhs) {
    auto counter = std::array<size_t, 2>{0, 0};

    return suffixarray::SparseArray<T>{
        R | std::views::transform([&](size_t b) {
            counter[b] += 1;
            if (!b) {
                return lhs.value(counter[0]-1);
            } else {
                return rhs.value(counter[1]-1);
            }
        })
    };
}
}

/**
 * Merges two FMIndices into a new one
 */
template <size_t SigmaLhs, template <size_t> typename StringLhs, size_t SigmaRhs, template <size_t> typename StringRhs, typename SuffixArray>
auto merge(FMIndex<SigmaLhs, StringLhs, SuffixArray> const& lhs, FMIndex<SigmaRhs, StringRhs, SuffixArray> const& rhs) -> FMIndex<std::max(SigmaLhs, SigmaRhs), StringLhs, SuffixArray> {
    auto R         = detail::computeInterleavingR(lhs.bwt, rhs.bwt);
    auto mergedBwt = detail::mergeBwt(R, lhs.bwt, rhs.bwt);
    auto annotatedArray = detail::mergeSparseArrays(R, lhs.annotatedArray, rhs.annotatedArray);

    return {mergedBwt, std::move(annotatedArray)};

}

/**
 * Merges two bidirectional FMIndices into a new one
 */
template <size_t SigmaLhs, template <size_t> typename StrLhs, size_t SigmaRhs, template <size_t> typename StrRhs, typename SuffixArray>
auto merge(BiFMIndex<SigmaLhs, StrLhs, SuffixArray> const& lhs, BiFMIndex<SigmaRhs, StrRhs, SuffixArray> const& rhs) -> BiFMIndex<std::max(SigmaLhs, SigmaRhs), StrLhs, SuffixArray> {
    // compute R of fwd BWT and merge bwt and csa
    auto R              = detail::computeInterleavingR(lhs.bwt, rhs.bwt);
    auto mergedBwt      = detail::mergeBwt(R, lhs.bwt, rhs.bwt);
    auto annotatedArray = detail::mergeSparseArrays(R, lhs.annotatedArray, rhs.annotatedArray);

    // compute R of rev BWT and merge BwtRev (reusing R to save on space)
    R                 = detail::computeInterleavingR(lhs.bwtRev, rhs.bwtRev);
    auto mergedBwtRev = detail::mergeBwt(R, lhs.bwtRev, rhs.bwtRev);

    return {mergedBwt, mergedBwtRev, std::move(annotatedArray)};
}
}
