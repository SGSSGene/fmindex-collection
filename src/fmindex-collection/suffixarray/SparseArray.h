// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../bitvector/L0L1_NBitvector.h"
#include "concepts.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <optional>
#include <tuple>

#if __has_include(<cereal/types/tuple.hpp>)
#include <cereal/types/tuple.hpp>
#endif
#if __has_include(<cereal/types/vector.hpp>)
#include <cereal/types/vector.hpp>
#endif


namespace fmindex_collection::suffixarray {

/** A sparse array
 *
 * Space efficient array over sparse data
 */
template <typename T>
struct SparseArray {
    using Bitvector = bitvector::L0L1_512_64kBitvector;
    std::vector<T> documentContent;
    Bitvector      bv;

    SparseArray() = default;
    SparseArray(SparseArray const&) = delete;
    SparseArray(SparseArray&&) noexcept = default;

    template <std::ranges::range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, std::optional<T>>
    SparseArray(range_t const& _range)
        : bv{_range | std::views::transform([&](std::optional<T> const& t) -> bool {
            if (t) {
                documentContent.push_back(*t);
                return true;
            }
            return false;
        })}
        {
    }

    auto operator=(SparseArray const&) -> SparseArray& = delete;
    auto operator=(SparseArray&&) noexcept -> SparseArray& = default;

    auto value(size_t idx) const -> std::optional<T> {
        assert(idx < bv.size());
        if (!bv.symbol(idx)) {
            return std::nullopt;
        }
        auto r = bv.rank(idx);
        assert(r < documentContent.size());
        auto v = documentContent[r];
        return v;
    }

    void push_back(std::optional<T> value) {
        bv.push_back(value.has_value());
        if (value) {
            documentContent.push_back(*value);
        }
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(documentContent, bv);
    }
};

}
