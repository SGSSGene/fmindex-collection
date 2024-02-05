// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <tuple>
#include <utility>
#include <vector>


namespace fmindex_collection {
/*
 * Minimum requirements to function as an Occurrence Table (OccTable)
 */
template<typename T, typename SymbolType = uint8_t>
concept OccTable = requires(T t, std::span<SymbolType const> bwt, uint64_t idx, SymbolType symb) {
    /** Every occtable has to be creatable by providing a bwt
     */
    { T{bwt} } -> std::same_as<T>;

    /** Every occtable has to have a C'Tor that can be constructed by default
     */
    { T{} } -> std::same_as<T>;

    /** Returns the name of this Occtable
     */
    { T::name() } -> std::same_as<std::string>;

    /** Returns the extension for this occtable
     */
    { T::extension() } -> std::same_as<std::string>;

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
    { t.all_ranks(idx) } -> std::same_as<std::tuple<std::array<uint64_t, T::Sigma>, std::array<uint64_t, T::Sigma>>>;

    /* Compile time variable indicating the number of symbols (including the delimiter)
     */
    { decltype(T::Sigma){} } -> std::same_as<size_t>;

    /* Run time variable indicating the number of rows inside this occurrence table
     */
    { t.size() } -> std::same_as<size_t>;

    /* Returns the symbol of the bwt at a certain position
     *
     * \param idx - row index
     * \return - will be in range of [0, Sigma)
     */
    //!TODO should not be equal to TLengthType
//    { t.symbol(idx) } -> (std::same_as<SymbolType> || std::same_as<TLengthType>);
};

template<template <auto> typename T>
concept checkOccTable = /*OccTable<T<1>>
                     && OccTable<T<2>>
                     &&*/ OccTable<T<4>>
                     && OccTable<T<5>>
                     && OccTable<T<254>>;
//                     && OccTable<T<256>>;



/** Additional prefetching methods that allows accelerated searches
 */
template <typename T>
concept OccTablePrefetch = OccTable<T> and requires(T t) {
    /*
     * Function that will prefetch the cache lines for
     * t.rank(...)
     * t.prefix_rank(...)
     * t.all_ranks(...)
     */
    t.prefetch(uint64_t{});
};


/** Additional methods to estimate memory usage
 */
template <typename T>
concept OccTableMemoryUsage = OccTable<T> and requires(T t) {
    /* Estimates the memory usage of a referce certain length
     *
     * \param first - length of the text to be indexed
     * \return number of bytes expected to be needed
     */
    { T::expectedMemoryUsage(uint64_t{}) } -> std::same_as<uint64_t>;

    /* Computes the actual memory usage of this data structure
     *
     * \return number of bytes that are being used
     */
    { t.memoryUsage() } -> std::same_as<uint64_t>;
};

}
