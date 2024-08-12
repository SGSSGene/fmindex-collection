#pragma once

#include <array>
#include <bitset>

namespace fmindex_collection {

template <size_t N>
inline std::array<std::bitset<N>, N+1> const leftshift_masks = []() {
    auto m = std::array<std::bitset<N>, N+1>{};
    m[0].flip();

    for (size_t i{1}; i < N+1; ++i) {
        m[i] = m[i-1] << 1;
    }
    return m;
}();

template <size_t N>
inline std::array<std::bitset<N>, N+1> const rightshift_masks = []() {
    auto m = std::array<std::bitset<N>, N+1>{};
    m[0].flip();

    for (size_t i{1}; i < N+1; ++i) {
        m[i] = m[i-1] >> 1;
    }
    return m;
}();

template <size_t N>
size_t lshift_and_count(std::bitset<N> const& b, size_t shift) {
    return (b & (leftshift_masks<N>[shift])).count();
}

template <size_t N>
size_t rshift_and_count(std::bitset<N> const& b, size_t shift) {
    return (b & (rightshift_masks<N>[shift])).count();
}


template <size_t N, typename Archive>
void loadBV(std::bitset<N>& b, Archive& ar) {
    b = std::bitset<N>{};
    for (size_t i{0}; i < N; i += 64) {
        uint64_t v{};
        ar(v);
        b = b | (std::bitset<N>{v} << i);
    }
}
template <size_t N, typename Archive>
void loadBV(std::vector<std::bitset<N>>& v, Archive& ar) {
    size_t ct{};
    ar(ct);
    v.resize(ct);
    for (size_t j{0}; j < ct; ++j) {
        loadBV(v[j], ar);
    }

}
template <size_t N, typename Archive>
void saveBV(std::bitset<N> const& b, Archive& ar) {
    (void)b;
    (void)ar;
    static constexpr auto mask = std::bitset<N>{~uint64_t{0}};

    // saving in 64bit blocks
    for (size_t i{0}; i < N; i += 64) {
        auto v = ((b >> i) & mask).to_ullong();
        ar(v);
    }
}
template <size_t N, typename Archive>
void saveBV(std::vector<std::bitset<N>> const& v, Archive& ar) {
    size_t ct = v.size();
    ar(ct);

    for (auto const& b : v) {
        saveBV(b, ar);
    }
}

}
