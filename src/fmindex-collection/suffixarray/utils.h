// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "SparseArray.h"

namespace fmindex_collection::suffixarray {

template <SuffixArray_c TCSA, SparseArray_c SparseArray = fmindex_collection::suffixarray::SparseArray<std::tuple<size_t, size_t>>>
auto convertCSAToAnnotatedDocument(TCSA const& _csa) -> SparseArray {
    return SparseArray{std::views::iota(size_t{0}, _csa.bv.size()) | std::views::transform([&](size_t i) {
        return _csa.value(i);
    })};
}

}
