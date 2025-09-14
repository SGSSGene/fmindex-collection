// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../ternarylogic.h"
#include "concepts.h"
#include "Bitvector.h"
#include "SparseBLEBitvector.h"
#include "L0L1_NBitvector.h"

#include <array>
#include <bitset>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ranges>
#include <stdexcept>
#include <span>
#include <vector>

#if __has_include(<cereal/types/variant.hpp>)
    #include <cereal/types/variant.hpp>
#endif

namespace fmc::bitvector {
/**
 * SparseDynRBBitvector sparse block length encoded bitvector.
 *
 * Uses Block length encoding and simplification on the
 * assumption the bit vector is sparse
 */
template <Bitvector_c BV1 = Bitvector, Bitvector_c BV2 = Bitvector>
struct SparseDynRBBitvector {

    using Variant = std::variant<
        L0L1_NBitvector<512ul, 65536ul>,
        SparseBLEBitvector<1, BV1, BV2>,
        SparseBLEBitvector<2, BV1, BV2>,
        SparseBLEBitvector<3, BV1, BV2>,
        SparseBLEBitvector<4, BV1, BV2>,
        SparseBLEBitvector<5, BV1, BV2>,
        SparseBLEBitvector<6, BV1, BV2>,
        SparseBLEBitvector<7, BV1, BV2>,
        SparseBLEBitvector<8, BV1, BV2>,
        SparseBLEBitvector<9, BV1, BV2>,
        SparseBLEBitvector<10, BV1, BV2>,
        SparseBLEBitvector<11, BV1, BV2>,
        SparseBLEBitvector<12, BV1, BV2>
    >;
    Variant bitvector;
    size_t totalLength{};


    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint8_t>
    SparseDynRBBitvector(range_t&& _range)
        : totalLength{_range.size()} {

        static constexpr size_t Level = std::variant_size_v<Variant>;
        auto countRuns = std::array<size_t, Level>{};

        uint64_t runCt{};
        for (size_t i{0}; i < _range.size(); ++i) {
            auto value = _range[i];
            runCt = (runCt+1) * (1-value);
            for (size_t j{1}; j < Level; ++j) {
                auto blockSize = (uint64_t{1} << j);
                if ((i % blockSize) == 0) {
                    if (runCt > blockSize) {
                        countRuns[j] += 1;
                    }
                }
            }
        }

        // initialized with the zeroth entry, which is a normal two layer 512bit bit vector
        size_t minElement = _range.size();
        size_t index{0};
        // Compute size of each:
        for (size_t j{1}; j < Level; ++j) {
            auto blockSize = (uint64_t{1} << j);
            auto indicatorEntries = _range.size() / blockSize;
            auto mixedBlockEntries = (indicatorEntries - countRuns[j]) * blockSize;
            auto totalSize = (indicatorEntries + mixedBlockEntries);
            if (totalSize < minElement) {
                minElement = totalSize;
                index = j;
            }
        }

        for_constexpr<0, Level>([&]<size_t L>() {
            if (L == index) {
                using V = std::variant_alternative_t<L, Variant>;
                bitvector = V{_range};
            }
        });
    }

    SparseDynRBBitvector() = default;
    SparseDynRBBitvector(SparseDynRBBitvector const&) = default;
    SparseDynRBBitvector(SparseDynRBBitvector&&) noexcept = default;

    auto operator=(SparseDynRBBitvector const&) -> SparseDynRBBitvector& = default;
    auto operator=(SparseDynRBBitvector&&) noexcept -> SparseDynRBBitvector& = default;

    size_t size() const noexcept {
        return totalLength;
    }

    bool symbol(size_t idx) const noexcept {
        return std::visit([&](auto const& bv) {
            return bv.symbol(idx);
        }, bitvector);
    }

    uint64_t rank(size_t idx) const noexcept {
        return std::visit([&](auto const& bv) {
            return bv.rank(idx);
        }, bitvector);
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bitvector, totalLength);
    }
};
static_assert(Bitvector_c<SparseDynRBBitvector<>>);

}
