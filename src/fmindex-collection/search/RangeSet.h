// SPDX-FileCopyrightText: 2026 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <map>
#include <cstdint>


namespace fmc {

/* A set of ranges
 */
class RangeSet {
public:
    std::map<size_t, size_t> ranges;
public:
    void clear() {
        ranges.clear();
    }

    // inserts the range starting at _pos with length _len
    void insert(size_t _pos, size_t _len) {
        auto iter = ranges.upper_bound(_pos+_len);
        // merge with all ranges that fall into this new range
        while (iter != ranges.begin()) {
            auto prev = std::prev(iter);
            auto const& [p, l] = *prev;
            if (p+l < _pos) break;
            auto newPos = std::min(p, _pos);
            auto newLen = std::max(p+l, _pos + _len) - newPos;
            _pos = newPos;
            _len = newLen;
            ranges.erase(prev);
            iter = ranges.upper_bound(_pos + _len);
        }

        ranges[_pos] = _len;
    }

    // checks if range [_pos, _len) is completely covered
    bool contains(size_t _pos, size_t _len) const {
        auto iter = ranges.upper_bound(_pos+_len);
        if (iter == ranges.begin()) return false;
        auto const& [pos, len] = *std::prev(iter);
        return _pos >= pos && _pos + _len <= pos+len;
    }

    // checks if point v is inside some range
    bool contains(size_t v) const {
        auto iter = ranges.upper_bound(v);
        if (iter == ranges.begin()) return false;
        auto const& [pos, len] = *std::prev(iter);
        return v < pos+len;
    }
};

}
