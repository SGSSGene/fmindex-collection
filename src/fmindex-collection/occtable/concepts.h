#pragma once

#include "../cereal_tag.h"

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
template<typename T, typename TLengthType = typename T::TLengthType>
concept OccTable =
    //!TODO stable on a single c'tor
    (
        requires(T t, std::vector<uint8_t> const& bwt, TLengthType idx, TLengthType symb) {
            { T{bwt} } -> std::same_as<T>;
        }
        || requires(T t, std::span<uint8_t const> bwt, TLengthType idx, TLengthType symb) {
            { T{bwt} } -> std::same_as<T>;
        }
    )
    && requires(T t, std::vector<uint8_t> const& bwt, TLengthType idx, TLengthType symb) {
    /** Every occtable has to be creatable by providing a bwt
     */
//    { T{bwt} } -> std::same_as<T>;

    /** Every occtable has to have a C'Tor that accept cereal_tag{}.
     * This constructor is used during deserialization
     */
    { T{cereal_tag{}} } -> std::same_as<T>;

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
    { t.rank(idx, symb) } -> std::same_as<TLengthType>;

    /* Return the numbers of symbols at a certain row that are equal or smaller.
     *
     * \param first - row index
     * \param second - symbol, a vale in the range of [1, Sigma)
     * \return number of occurrences
     */
    { t.prefix_rank(idx, symb) } -> std::same_as<TLengthType>;

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
    { t.all_ranks(idx) } -> std::same_as<std::tuple<std::array<TLengthType, T::Sigma>, std::array<TLengthType, T::Sigma>>>;

    /* Compile time variable indicating the number of symbols (including the delimiter)
     */
    { decltype(T::Sigma){} } -> std::same_as<TLengthType>;

    /* Run time variable indicating the number of rows inside this occurrence table
     */
    { t.size() } -> std::same_as<TLengthType>;

    /* Returns the symbol of the bwt at a certain position
     *
     * \param idx must be in range of [0, Sigma)
     */
    { t.symbol(idx) } -> std::same_as<TLengthType>;
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
