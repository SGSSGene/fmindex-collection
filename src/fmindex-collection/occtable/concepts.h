#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <tuple>
#include <utility>

template<typename T>
concept OccTable = requires(T t) {
    { t.rank(uint64_t{}, uint64_t{}) } -> uint64_t;
    { t.prefix_rank(uint64_t{}, uint64_t{}) } -> uint64_t;
    { t.all_ranks(uint64_t{}) } -> std::tuple<std::array<uint64_t, T::Sigma>, std::array<uint64_t, T::Sigma>>;
    { T::Sigma } -> uint64_t;

    t.size();
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


template <typename T>
concept OccTablePrefetch = OccTable<T> and requires(T t) {
    t.prefetch(uint64_t{});
};
