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

}
