// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../search_scheme/generator/h2.h"
#include "../search_scheme/expand.h"

namespace fmc {

/** cache a search scheme (without expansion to length)
 *
 * @param _shortLen: indicates that the search scheme should work for short queries (length of 2)
 */
template <bool Edit>
auto getCachedSearchScheme(size_t _minError, size_t _maxError, bool _shortLen=false) -> auto const& {
    static thread_local auto cache = std::tuple<size_t, size_t, bool, fmc::search_scheme::Scheme>{std::numeric_limits<size_t>::max(), 0, false, {}};
    // check if last scheme has correct errors and length, otherwise generate it
    auto& [minError, maxError, shortLen, search_scheme] = cache;
    if (_minError != minError
        || _maxError != maxError
        || shortLen != _shortLen
    ) { // regenerate everything
        minError      = _minError;
        maxError      = _maxError;
        shortLen      = _shortLen;
        search_scheme = fmc::search_scheme::generator::h2(maxError+(shortLen?1:2), minError, maxError);

        if constexpr (!Edit) {
            search_scheme = limitToHamming(search_scheme);
        }
    }
    return search_scheme;
}

/** cache a search scheme with expansion to length
 */
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

inline auto getCachedPartition(size_t _parts, size_t _length) -> auto const & {
    static thread_local auto cache = std::tuple<size_t, size_t, std::vector<size_t>>{0, 0, {}};
    auto& [length, parts, partition] = cache;
    if (length != _length || parts != _parts) {
        partition = fmc::search_scheme::createUniformPartition(_parts, _length);
        parts = _parts;
        length = _length;
    }
    return partition;
}

}
