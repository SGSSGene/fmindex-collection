// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../bitvector/concepts.h"
#include "../string/concepts.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <tuple>
#include <utility>

#if __has_include(<cereal/types/vector.hpp>)
#include <cereal/archives/binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/vector.hpp>
#endif


namespace fmc {

template <typename T, typename SymbolType = uint8_t>
concept String_c =
    std::default_initializable<T>
    && std::movable<T>
    && requires(T const t, std::span<SymbolType const> symbols, size_t idx, SymbolType symb) {

    /* Compile time variable indicating the number of symbols (including the delimiter)
     */
    { decltype(T::Sigma){} } -> std::same_as<size_t>;

    /** Every String_c can be constructed via some type of string similar thing
     */
    { T{symbols} } -> std::same_as<T>;

    /* Run time variable indicating the number of rows inside this occurrence table
     */
    { t.size() } -> std::same_as<size_t>;

    /* Returns the symbol of the symbol at a certain position
     *
     * \param idx - row index
     * \return - will be in range of [0, Sigma)
     */
    { t.symbol(idx) };// -> std::same_as<SymbolType>; //!TODO

    /* Return the numbers of symbols at a certain row + all symbols smaller over all rows
     *
     * \param first - row index
     * \param second - symbol, a value in the range of [1, Sigma)
     * \return number of occurrences
     */
    { t.rank(idx, symb) } -> std::same_as<uint64_t>;

    /* Return the numbers of symbols at a certain row that are equal or smaller.
     *
     * \param first - row index
     * \param second - symbol, a vale in the range of [1, Sigma)
     * \return number of occurrences
     */
    { t.prefix_rank(idx, symb) } -> std::same_as<uint64_t>;

    /* Combined rank and prefix_rank over all symbols
     * !TODO documentation outdated
     *
     * \param first - row index
     * \return a tuple of two arrays.
     *         The first array are the ranks (as returned by t.rank)
     *         The second array are the prefix ranks (as returned by t.prefix_ranks)
     * \example
     * auto [ranks, prefixes] = index.all_ranks(10);
     * assert(ranks[1] == index.rank(10, 1);
     * assert(ranks[2] == index.rank(10, 2);
     * \\ ...
     * assert(prefixes[1] == index.prefix_rank(10, 1);
     * assert(prefixes[1] == index.prefix_rank(10, 2);
     * \\ ...
     */
    { t.all_ranks(idx) } -> std::same_as<std::array<uint64_t, T::Sigma>>;

    /* !TODO needs documentation
     */
    { t.all_ranks_and_prefix_ranks(idx) } -> std::same_as<std::tuple<std::array<uint64_t, T::Sigma>, std::array<uint64_t, T::Sigma>>>;
};


template<template <auto> typename T>
concept checkString_c =
    String_c<T<2>>
    && String_c<T<4>>
    && String_c<T<5>>
    && String_c<T<255>>
    && String_c<T<256>>
;

template <typename T, typename SymbolType = uint8_t>
concept StringKStep_c = String_c<T, SymbolType>
    && requires(T const t, std::span<SymbolType const> symbols, size_t idx, SymbolType symb) {
    /* Returns the symbol at a certain position only looking at the last L bits
     *
     * \tparam L  - number of bits at the end
     * \param idx - row index
     * \return - will be in range of [0, (1<<L))
     */
    { t.symbol_limit(idx) };// -> std::same_as<SymbolType>; //!TODO

    /* Same as rank, but only looks at the last L bits
     */
    { t.rank_limit(idx, symb) } -> std::same_as<uint64_t>;

    /* Fetches prefix rank and rank for the last L bits
     */
    { t.prefix_rank_and_rank_limit(idx, symb) } -> std::same_as<std::tuple<uint64_t, uint64_t>>;

};

template<template <auto> typename T>
concept checkStringKStep_c =
    StringKStep_c<T<2>>
    && StringKStep_c<T<4>>
    && StringKStep_c<T<5>>
    && StringKStep_c<T<255>>
    && StringKStep_c<T<256>>
;


}
