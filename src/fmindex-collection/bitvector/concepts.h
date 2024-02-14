// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <tuple>
#include <utility>
#include <type_traits>

#if __has_include(<cereal/types/vector.hpp>)
#include <cereal/archives/binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/vector.hpp>
#endif

namespace fmindex_collection {

template <typename T>
concept BitVector_c =
    std::default_initializable<T>
    && std::movable<T>
    && requires(T t, std::span<uint8_t const> symbols, size_t idx, uint8_t symb) {

    /** Every BitVector can be constructed via some type of string similar thing
     */
    { T{size_t{}, [](size_t) { return false; }} } -> std::same_as<T>;

    /* Run time variable indicating the number of rows inside this occurrence table
     */
    { t.size() } -> std::same_as<size_t>;

    /* Returns the symbol of the symbol at a certain position
     *
     * \param idx - row index
     * \return - will be true or false
     */
    { t.symbol(idx) } -> std::same_as<bool>;

    /* Return the numbers of ones
     *
     * \param first - row index
     * \return number of occurrences
     */
    { t.rank(idx) } -> std::same_as<uint64_t>;
};

}
