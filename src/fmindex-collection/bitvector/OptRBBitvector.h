// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../ternarylogic.h"
#include "concepts.h"
#include "Bitvector.h"
#include "RBBitvector.h"
#include "Bitvector2L.h"

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
 * OptRBBitvector sparse block length encoded bitvector.
 *
 * Uses Block length encoding and simplification on the
 * assumption the bit vector is sparse
 */
template <Bitvector_c BV1 = Bitvector, Bitvector_c BV2 = Bitvector>
struct OptRBBitvector {

    using Variant = std::variant<
        Bitvector2L<512ul, 65536ul>,
        RBBitvector<1, BV1, BV2>,
        RBBitvector<2, BV1, BV2>,
        RBBitvector<3, BV1, BV2>,
        RBBitvector<4, BV1, BV2>,
        RBBitvector<5, BV1, BV2>,
        RBBitvector<6, BV1, BV2>,
        RBBitvector<7, BV1, BV2>,
        RBBitvector<8, BV1, BV2>,
        RBBitvector<9, BV1, BV2>,
        RBBitvector<10, BV1, BV2>,
        RBBitvector<11, BV1, BV2>,
        RBBitvector<12, BV1, BV2>
    >;

    Variant bitvector;
    size_t totalLength{};


    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint8_t>
    OptRBBitvector(range_t&& _range)
        : totalLength{_range.size()} {

        static constexpr size_t Level = std::variant_size_v<Variant>;
        auto countRuns = std::array<std::tuple<size_t, size_t>, Level>{};

        uint64_t runCtZero{};
        uint64_t runCtOne{};
        for (size_t i{0}; i < _range.size(); ++i) {
            auto value = _range[i];
            runCtZero = (runCtZero+1) * (1-value);
            runCtOne = (runCtOne+1) * value;

            for_constexpr<1, Level>([&]<size_t L>() {
                using V = std::variant_alternative_t<L, Variant>;
                if ((i % V::BlockLength) != 0) return;
                if (runCtZero > V::BlockLength) {
                    std::get<0>(countRuns[L]) += 1;
                }
                if (runCtOne > V::BlockLength) {
                    std::get<1>(countRuns[L]) += 1;
                }
            });
        }

        // initialized with the zeroth entry, which is a normal two layer 512bit bit vector
        size_t minElement = std::variant_alternative_t<0, Variant>::estimateSize(_range.size());
//        size_t minElement = std::numeric_limits<size_t>::max();
        size_t index{0};
        // Compute size of each:
        for_constexpr<1, Level>([&]<size_t L>() {
            using V = std::variant_alternative_t<L, Variant>;
            auto totalSize = V::estimateSize(_range.size() / V::BlockLength, std::get<0>(countRuns[L]), std::get<1>(countRuns[L]));
            if (totalSize < minElement) {
                minElement = totalSize;
                index = L;
            }
        });

        for_constexpr<0, Level>([&]<size_t L>() {
            if (L == index) {
                using V = std::variant_alternative_t<L, Variant>;
                bitvector = V{_range};
            }
        });
    }

    OptRBBitvector() = default;
    OptRBBitvector(OptRBBitvector const&) = default;
    OptRBBitvector(OptRBBitvector&&) noexcept = default;

    auto operator=(OptRBBitvector const&) -> OptRBBitvector& = default;
    auto operator=(OptRBBitvector&&) noexcept -> OptRBBitvector& = default;

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
static_assert(Bitvector_c<OptRBBitvector<>>);

}
