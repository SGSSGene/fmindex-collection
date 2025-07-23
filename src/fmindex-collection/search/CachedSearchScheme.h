// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../search_scheme/generator/h2.h"
#include "../search_scheme/expand.h"

namespace fmc {

template <bool Edit>
auto getCachedSearchScheme(size_t _length, size_t _minError, size_t _maxError) -> auto const& {
    static thread_local auto cache = std::tuple<size_t, size_t, size_t, fmc::search_scheme::Scheme, fmc::search_scheme::Scheme>{std::numeric_limits<size_t>::max(), 0, 0, {}, {}};
    // check if last scheme has correct errors and length, other wise generate it
    auto& [minError, maxError, length, ss, search_scheme] = cache;
    if (_minError != minError || _maxError != maxError) { // regenerate everything
        minError      = _minError;
        maxError      = _maxError;
        ss            = fmc::search_scheme::generator::h2(maxError+2, minError, maxError);
        length        = _length;
        search_scheme = fmc::search_scheme::expand(ss, length);
        if constexpr (!Edit) {
            search_scheme = limitToHamming(search_scheme);
        }
    } else if (_length != length) {
        length        = _length;
        search_scheme = fmc::search_scheme::expand(ss, length);
        if constexpr (!Edit) {
            search_scheme = limitToHamming(search_scheme);
        }
    }
    return search_scheme;
}

}
