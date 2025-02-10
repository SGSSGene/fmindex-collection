// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace fmindex_collection {

/* A dense std::vector<uint> implementation using a fixed width (number of bits).
 *
 * This uses a fixed width/a fixed number of bits for every entry. Allowing for a dense implementation of
 * vectors over integers.
 *
 * This vector only implements a subset of `std::vector` functionality.
 *
 */
template <size_t WordWidth, size_t MaxValue = std::numeric_limits<size_t>::max()>
struct FixedSuccinctVector {
    std::vector<uint64_t> data; // buffer where the data is being stored
    size_t bitCount{};          // numbers of used bits

    FixedSuccinctVector() = default;
    FixedSuccinctVector(std::initializer_list<uint64_t> list) {
        for (auto s : list) {
            push_back(s);
        }
    }
    FixedSuccinctVector(FixedSuccinctVector const&) = default;
    FixedSuccinctVector(FixedSuccinctVector&&) noexcept = default;
    FixedSuccinctVector& operator=(FixedSuccinctVector const&) = default;
    FixedSuccinctVector& operator=(FixedSuccinctVector&&) noexcept = default;

    /** Reserve memory for a certain number of integers.
     *
     * \param s number of integer
     */
    void reserve(size_t s) {
        auto a = WordWidth*s;

        data.reserve(a/64+((a % 64)>0?1:0));
    }


    /** Push a integer at the end of the vector
     *
     * \param value value to push onto the vector. Only the first `WordWidth` are
     * being stored.
     */
    void push_back(uint64_t value) {
        assert(value <= MaxValue);
        assert(std::log2(value) < WordWidth);
        auto empty = data.size()*64-bitCount;
        if (empty == 0) {
            data.push_back(value);
            bitCount += WordWidth;
            assert(back() == value);
            return;
        }
        if (empty >= WordWidth) {
            data.back() |= value << (64-empty);
            bitCount += WordWidth;
            assert(back() == value);
            return;
        }

        data.back() |= value << (64-empty); // only pushes some bits, some are being dropped
        data.push_back(value >> empty);
        bitCount += WordWidth;
        assert(back() == value);
    }
    void emplace_back(uint64_t value) {
        push_back(value);
    }

    /** Read integer at a certain position
     */
    auto operator[](size_t i) const {
        return at(i);
    }

    /** Read integer at a certain position
     */
    auto at(size_t i) const -> uint64_t {
        auto begin = i * WordWidth;
        auto end   = begin+WordWidth-1; // The end is inclusive
        assert(begin == end || begin < bitCount);
        assert(end <= bitCount);
        assert(end >= begin);
        assert(end-begin <= 64);

        auto startI = begin/64;
        auto endI   = end/64;
        auto startOffset = begin % 64;

        auto mask = [&]() -> uint64_t {
            if constexpr (WordWidth < 64) {
                return (uint64_t{1}<<(end-begin+size_t{1}))-1;
            } else if constexpr (WordWidth == 64) {
                return std::numeric_limits<uint64_t>::max();
            }
        }();

        if (startI == endI) {
            return (data[startI] >> startOffset) & mask;
        }

        auto p1 = data[startI] >> startOffset;
        auto p2 = data[endI] << (64-startOffset);
        auto v = (p1 | p2) & mask;
        assert(std::log2(v) < WordWidth);
        return v;
    }

    void pop_back() {
        bitCount -= WordWidth;
        if (data.size()*64 - bitCount >= 64) {
            data.pop_back();
        }
        // clear left over zeros
        auto empty = data.size()*64-bitCount;
        if (empty == 0) return;
        auto mask = std::numeric_limits<uint64_t>::max() >> empty;
        data.back() = data.back() & mask;
    }

    uint64_t back() const {
        return at(size()-1);
    }

    size_t size() const {
        return bitCount / WordWidth; // Always a whole number
    }

    /** Number of bits required to store all integers
     */
    size_t bit_size() const {
        return bitCount;
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(data, bitCount);
    }
};

}
