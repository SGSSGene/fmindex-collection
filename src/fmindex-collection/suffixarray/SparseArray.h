// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../bitvector/Bitvector2L.h"
#include "../bitvector/PairedBitvector2L.h"
#include "concepts.h"
#include "../DenseVector.h"

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

template <typename T>
struct SparseArray;


/** A sparse array over multiple values
 *
 * Space efficient array over sparse data
 */
template <typename... Ts>
    requires (std::convertible_to<Ts, size_t> && ...)
        && (std::convertible_to<size_t, Ts> && ...)
struct SparseArray<std::tuple<Ts...>> {
    using Entry     = std::tuple<Ts...>;
    using value_t   = Entry;
    using Bitvector = bitvector::Bitvector2L<512, 65536>;
    using Documents = std::array<DenseVector, sizeof...(Ts)>;
    Documents documents;
    Bitvector bv;

    SparseArray() = default;
    SparseArray(SparseArray const&) = delete;
    SparseArray(SparseArray&&) noexcept = default;

    template <size_t End>
    static void map(auto const& cb) {
        if constexpr (End > 0) {
            cb.template operator()<End-1>();
            map<End-1>(cb);
        }
    }

    template <std::ranges::range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, std::optional<Entry>>
    SparseArray(range_t const& _range) {
        auto largestValue = std::array<size_t, sizeof...(Ts)>{};
        auto commonDivisor = std::array<size_t, sizeof...(Ts)>{};
        size_t ct{};

        bv = Bitvector{_range | std::views::transform([&](std::optional<Entry> const& t) -> bool {
            if (t) {
                // using template magic to spread a single entry into multiple vectors
                map<sizeof...(Ts)>([&]<size_t I>() {
                    auto const v = static_cast<size_t>(std::get<I>(*t));
                    largestValue[I]  = std::max(largestValue[I], v);
                    commonDivisor[I] = std::gcd(commonDivisor[I], v);
                });
                ct += 1;
                return true;
            }
            return false;
        })};
        map<sizeof...(Ts)>([&]<size_t I>() {
            if (commonDivisor[I] == 0) {
                commonDivisor[I] = 1;
            }
            if (largestValue[I] == 0) {
                largestValue[I] = 1;
            }
        });


        map<sizeof...(Ts)>([&]<size_t I>() {
            documents[I] = DenseVector(/*.largestValue = */ largestValue[I], /*.commonDivisor = */commonDivisor[I]);
            documents[I].reserve(ct);
        });
        for (auto t : _range) {
            if (!t) continue;
            map<sizeof...(Ts)>([&]<size_t I>() {
                auto const v = static_cast<size_t>(std::get<I>(*t));
                documents[I].push_back(v);
            });
        }
    }

    auto operator=(SparseArray const&) -> SparseArray& = delete;
    auto operator=(SparseArray&&) noexcept -> SparseArray& = default;

    auto value(size_t idx) const -> std::optional<Entry> {
        assert(idx < bv.size());
        if (!bv.symbol(idx)) {
            return std::nullopt;
        }
        auto r = bv.rank(idx);
        auto entry = Entry{};
        // using template magic to gather multiple vectors into a single entry
        map<sizeof...(Ts)>([&]<size_t I>() {
            assert(r < std::get<I>(documents).size());
            std::get<I>(entry) = std::get<I>(documents)[r];
        });

        return entry;
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(documents, bv);
    }
};

}
