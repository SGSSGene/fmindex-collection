// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <cstdint>
#include <ranges>

namespace fmindex_collection {

template <typename T>
concept Sequence  = std::ranges::sized_range<T>
                    && std::ranges::random_access_range<T>
                    && requires(T t) {
                        {*t.begin()} -> std::common_with<uint8_t>;
                    };


template <typename T>
concept Sequences = std::ranges::sized_range<T>
                    && std::ranges::random_access_range<T>
                    && requires(T t) {
                        {*t.begin()} -> Sequence;
                    };

auto constexpr add_sentinel = std::views::transform([] (auto r) -> uint8_t {
    return r + 1;
});

auto constexpr add_sentinels = std::views::transform([] (auto seq) {
    return seq | add_sentinel;
});

template <typename T>
concept Sequence_32  = std::ranges::sized_range<T>
                    && std::ranges::random_access_range<T>
                    && requires(T t) {
                        {*t.begin()} -> std::common_with<uint32_t>;
                    };


template <typename T>
concept Sequences_32 = std::ranges::sized_range<T>
                    && std::ranges::random_access_range<T>
                    && requires(T t) {
                        {*t.begin()} -> Sequence_32;
                    };


}
