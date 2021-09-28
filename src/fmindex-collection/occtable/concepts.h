#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <tuple>
#include <utility>

/*
 * Minimum requirements to function as an Occurrence Table (OccTable)
 */
template<typename T>
concept OccTable = requires(T t) {
    /* Return the numbers of symbols at a certain row + all symbols smaller over all rows
     *
     * \param first - row index
     * \param second - symbol, a value in the range of [1, Sigma)
     * \return number of occurrences
     */
    { t.rank(uint64_t{}, uint64_t{}) } -> uint64_t;

    /* Return the numbers of symbols at a certain row that are equal or smaller.
     *
     * \param first - row index
     * \param second - symbol, a vale in the range of [1, Sigma)
     * \return number of occurrences
     */
    { t.prefix_rank(uint64_t{}, uint64_t{}) } -> uint64_t;

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
    { t.all_ranks(uint64_t{}) } -> std::tuple<std::array<uint64_t, T::Sigma>, std::array<uint64_t, T::Sigma>>;

    /* Compile time variable indicating the number of symbols (including the delimiter)
     */
    { T::Sigma } -> uint64_t;

    /* Run time variable indicating the number of rows inside this occurrence table
     */
    { t.size() } -> size_t;

    T::expectedMemoryUsage(0);
    t.symbol(uint64_t{});
    t.memoryUsage();
};

template<template <auto> typename T>
concept checkOccTable = OccTable<T<1>>
                     && OccTable<T<2>>
                     && OccTable<T<4>>
                     && OccTable<T<5>>
                     && OccTable<T<254>>;



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
