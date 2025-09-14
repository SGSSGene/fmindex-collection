// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../bitvector/Bitvector2L.h"
#include "../bitvector/PairedBitvector2L.h"
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


namespace fmc::suffixarray {

/** A sparse array
 *
 * Space efficient array over sparse data
 */
template <typename T>
struct SparseArray {
    using value_t = T;
    using Bitvector = bitvector::Bitvector2L<512, 65536>;
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

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(documentContent, bv);
    }
};

}
