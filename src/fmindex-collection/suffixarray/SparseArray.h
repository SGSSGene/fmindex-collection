// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../bitvector/Bitvector2L.h"
#include "../bitvector/PairedBitvector2L.h"
#include "../bitvector/OptSparseRBBitvector.h"
#include "concepts.h"
#include "../DenseMultiVector.h"

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

/** A sparse array over multiple values
 *
 * Space efficient array over sparse data
 */
template <typename TEntry, typename TBitvector = bitvector::Bitvector2L<512, 65536>>
struct SparseArray {
    using Entry     = TEntry;
    using value_t   = Entry;
    using Bitvector = TBitvector;
    using Documents = DenseMultiVector<Entry>;
    Documents documents;
    Bitvector bv;

    SparseArray() = default;
    SparseArray(SparseArray const&) = delete;
    SparseArray(SparseArray&&) noexcept = default;

    template <std::ranges::range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, std::optional<Entry>>
    SparseArray(range_t const& _range) {
        bv = Bitvector{_range | std::views::transform([](std::optional<Entry> const& t) -> bool {
            return t.has_value();
        })};


        auto data = std::vector<Entry>{};
        for (auto const& v : _range) {
            if (!v) continue;
            data.push_back(*v);
        }
        documents = Documents { data };
    }

    auto operator=(SparseArray const&) -> SparseArray& = delete;
    auto operator=(SparseArray&&) noexcept -> SparseArray& = default;

    auto value(size_t idx) const -> std::optional<Entry> {
        assert(idx < bv.size());
        if (!bv.symbol(idx)) {
            return std::nullopt;
        }
        auto r = bv.rank(idx);
        return documents[r];
    }

    template <typename Archive>
    void serialize(this auto&& self, Archive& ar) {
        ar(self.documents, self.bv);
    }
};

}
