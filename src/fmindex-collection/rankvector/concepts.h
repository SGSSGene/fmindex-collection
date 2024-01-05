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

namespace fmindex_collection {

template <typename T>
concept BitVector = requires(T t, std::span<uint8_t const> symbols, size_t idx, uint8_t symb) {
    /** Every BitVector can be constructed via some type of string similar thing
     */
    { T{symbols} } -> std::same_as<T>;

    /** Every BitVector can be constructed via some type of string similar thing
     */
    { T{symbols, symb} } -> std::same_as<T>;

    /** Every RankVector has a C'Tor that accept cereal_tag{}.
     * This constructor is used during deserialization
     **/
    { T{cereal_tag{}} } -> std::same_as<T>;

    /* Run time variable indicating the number of rows inside this occurrence table
     */
    { t.size() } -> std::same_as<size_t>;

    /* Returns the symbol of the symbol at a certain position
     *
     * \param idx - row index
     * \return - will be in range of [0, Sigma)
     */
    { t.symbol(idx) } -> std::same_as<bool>;

    /* Return the numbers of ones
     *
     * \param first - row index
     * \param second - symbol, a value in the range of [1, Sigma)
     * \return number of occurrences
     */
    { t.rank(idx) } -> std::same_as<size_t>;

    /* Compile time variable indicating the number of symbols (including the delimiter)
     */
    { decltype(T::Sigma){} } -> std::same_as<size_t>;
};


template<typename T>
concept checkBitVector = BitVector<T>;

template <typename T, typename SymbolType = uint8_t>
concept SymbolVector = requires(T t, std::span<SymbolType const> symbols, size_t idx, SymbolType symb) {
    /* Compile time variable indicating the number of symbols (including the delimiter)
     */
    { decltype(T::Sigma){} } -> std::same_as<size_t>;

    /** Every RankVector can be constructed via some type of string similar thing
     */
    { T{symbols} } -> std::same_as<T>;

    /** Every RankVector has a C'Tor that accept cereal_tag{}.
     * This constructor is used during deserialization
     **/
    { T{cereal_tag{}} } -> std::same_as<T>;

    /* Run time variable indicating the number of rows inside this occurrence table
     */
    { t.size() } -> std::same_as<size_t>;

    /* Returns the symbol of the symbol at a certain position
     *
     * \param idx - row index
     * \return - will be in range of [0, Sigma)
     */
    { t.symbol(idx) } -> std::same_as<SymbolType>;

    /* Return the numbers of symbols at a certain row + all symbols smaller over all rows
     *
     * \param first - row index
     * \param second - symbol, a value in the range of [1, Sigma)
     * \return number of occurrences
     */
    { t.rank(idx, symb) } -> std::same_as<size_t>;

    /* Return the numbers of symbols at a certain row that are equal or smaller.
     *
     * \param first - row index
     * \param second - symbol, a vale in the range of [1, Sigma)
     * \return number of occurrences
     */
    { t.prefix_rank(idx, symb) } -> std::same_as<size_t>;

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
    { t.all_ranks(idx) } -> std::same_as<std::array<size_t, T::Sigma>>;

    /* !TODO needs documentation
     */
    { t.all_ranks_and_prefix_ranks(idx) } -> std::same_as<std::tuple<std::array<size_t, T::Sigma>, std::array<size_t, T::Sigma>>>;
};


template<template <auto> typename T>
concept checkSymbolVector = SymbolVector<T<1>>
                         && SymbolVector<T<2>>
                         && SymbolVector<T<4>>
                         && SymbolVector<T<5>>
                         && SymbolVector<T<255>>
                         && SymbolVector<T<256>>;
}
