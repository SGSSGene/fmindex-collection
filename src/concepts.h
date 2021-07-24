#pragma once

#include <utility>
#include <cstdint>
#include <cstddef>

template<typename T>
concept FMIndex = requires(T t) {
    t.rank(uint8_t{}, uint64_t{});
    t.prefix_rank(uint8_t{}, uint64_t{});
    T::Sigma;
    t.size();
    T::expectedMemoryUsage(0);
};

template<template <auto> typename T>
concept checkFMIndex = FMIndex<T<1>>
                    && FMIndex<T<2>>
                    && FMIndex<T<4>>
                    && FMIndex<T<5>>
                    && FMIndex<T<254>>;
