// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "DenseVector.h"

#include <array>
#include <variant>

namespace fmc {

template <typename Ts>
struct DenseMultiVector;

template <>
struct DenseMultiVector<std::monostate> {
    DenseMultiVector() = default;

    static auto concat(DenseMultiVector const& lhs, DenseMultiVector const& rhs) {
        (void)lhs;
        (void)rhs;
        return DenseMultiVector{};
    }

    DenseMultiVector(DenseMultiVector const&) = default;
    DenseMultiVector(DenseMultiVector&&) noexcept = default;
    DenseMultiVector& operator=(DenseMultiVector const&) = default;
    DenseMultiVector& operator=(DenseMultiVector&&) noexcept = default;

    void reserve(size_t) {}
    void push_back(std::monostate) {}
    auto operator[](size_t) const -> std::monostate { return {}; }
    auto access(size_t) const -> std::monostate { return {}; }
    size_t size() const { return 0; }

    template <typename Archive>
    void serialize(this auto&& self, Archive& ar) {
        (void)self;
        (void)ar;
    }
};


/* Like DenseVector but stores multiple values per entry
 */
template <typename ...Ts>
    requires (std::convertible_to<Ts, size_t> && ...)
        && (std::convertible_to<size_t, Ts> && ...)
struct DenseMultiVector<std::tuple<Ts...>> {
    std::array<DenseVector, sizeof...(Ts)> data;

    DenseMultiVector() = default;

    static auto concat(DenseMultiVector const& lhs, DenseMultiVector const& rhs) {
        auto v = DenseMultiVector{};
        for (size_t i{0}; i < lhs.data.size(); ++i) {
            v.data[i] = DenseVector::concat(lhs.data[i], rhs.data[i]);
        }
        return v;
    }


    DenseMultiVector(DenseMultiVector const&) = default;
    DenseMultiVector(DenseMultiVector&&) noexcept = default;
    DenseMultiVector& operator=(DenseMultiVector const&) = default;
    DenseMultiVector& operator=(DenseMultiVector&&) noexcept = default;

    template <std::ranges::range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, std::tuple<Ts...>>
    DenseMultiVector(range_t const& _range) {
        auto largestValue = std::array<size_t, sizeof...(Ts)>{};
        auto commonDivisor = std::array<size_t, sizeof...(Ts)>{};
        size_t ct{};

        // find largest value
        for (auto const& t : _range) {
            // using template magic to spread a single entry into multiple vectors
            map<sizeof...(Ts)>([&]<size_t I>() {
                auto const v = static_cast<size_t>(std::get<I>(t));
                largestValue[I]  = std::max<size_t>(largestValue[I], v);
                commonDivisor[I] = std::gcd(commonDivisor[I], v);
            });
            ct += 1;
        }

        // make sure commonDivisor and largestValue are at least 1
        map<sizeof...(Ts)>([&]<size_t I>() {
            if (commonDivisor[I] == 0) {
                commonDivisor[I] = 1;
            }
            if (largestValue[I] == 0) {
                largestValue[I] = 1;
            }
        });

        // Initialize each vector accordinlgy
        map<sizeof...(Ts)>([&]<size_t I>() {
            data[I] = DenseVector(/*.largestValue = */ largestValue[I], /*.commonDivisor = */commonDivisor[I]);
            data[I].reserve(ct);
        });
        for (auto const& t : _range) {
            push_back(t);
        }
    }


    /** Reserve memory for a certain number of integers.
     *
     * \param s number of integer
     */
    void reserve(size_t s) {
        for (size_t i{0}; i < data.size(); ++i) {
            data[i].reserve(s);
        }
    }

    template <size_t End>
    static void map(auto const& cb) {
        if constexpr (End > 0) {
            cb.template operator()<End-1>();
            map<End-1>(cb);
        }
    }

    /** Push a value to the end of the vector
     */
    void push_back(std::tuple<Ts...> const& value) {
        map<sizeof...(Ts)>([&]<size_t I>() {
            data[I].push_back(std::get<I>(value));
        });
    }

    /** Read integer at a certain position
     */
    auto operator[](size_t i) const -> std::tuple<Ts...> {
        return access(i);
    }

    /** Read integer at a certain position
     */
    auto access(size_t i) const -> std::tuple<Ts...> {
        auto ret = std::tuple<Ts...>{};
        map<sizeof...(Ts)>([&]<size_t I>() {
            std::get<I>(ret) = data[I].access(i);
        });
        return ret;
    }

    size_t size() const {
        return data[0].size();
    }

    template <typename Archive>
    void serialize(this auto&& self, Archive& ar) {
        ar(self.data);
    }
};

}
