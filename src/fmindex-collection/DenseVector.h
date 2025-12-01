// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <bit>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <span>
#include <vector>

namespace fmc {

/* A dense std::vector<uint> implementation using a fixed width (number of bits).
 *
 * This uses a fixed width/a fixed number of bits for every entry. Allowing for a dense implementation of
 * vectors over integers.
 *
 * This vector only implements a subset of `std::vector` functionality.
 *
 */
struct DenseVector {
    std::vector<uint64_t> data; // buffer where the data is being stored
    size_t  bitCount{};          // numbers of used bits
    uint8_t bits{};              // number of bits per entry
    size_t  largestValue{};      // largestValue that has to be stored
    size_t  commonDivisor{1};    // factor each value is multiplied with

    DenseVector() = default;


    static auto concat(DenseVector const& lhs, DenseVector const& rhs) {
        auto v = DenseVector{};
        v.largestValue  = std::max(lhs.largestValue, rhs.largestValue);
        v.commonDivisor = std::gcd(lhs.commonDivisor, rhs.commonDivisor);
        v.bits          = std::bit_width(v.largestValue / v.commonDivisor);
        for (size_t i{0}; i < lhs.size(); ++i) {
            v.push_back(lhs[i]);
        }
        for (size_t i{0}; i < rhs.size(); ++i) {
            v.push_back(rhs[i]);
        }
        return v;
    }

    /**
     * \param _bits number of bits used for each integer entry
     */
    DenseVector(size_t _largestValue, size_t _commonDivisor = 1)
        : bits{uint8_t(std::bit_width(_largestValue / _commonDivisor))}
        , largestValue{_largestValue}
        , commonDivisor{_commonDivisor}
    {}

    /**
     * Create a DenseVector and auto determine the number of bits required
     */
    DenseVector(std::initializer_list<uint64_t const> args)
        : commonDivisor{0}
    {
        for (auto v : args) {
            largestValue = std::max(largestValue, v);
            commonDivisor = std::gcd(commonDivisor, v);
        }
        if (commonDivisor == 0) {
            commonDivisor = 1;
        }

        bits = std::bit_width(largestValue / commonDivisor);
        for (auto v : args) {
            push_back(v);
        }
    }

    /**
     * Create a DenseVector and auto determine the number of bits required
     */
    DenseVector(std::span<uint64_t const> args)
        : commonDivisor{0}
    {
        for (auto v : args) {
            largestValue = std::max(largestValue, v);
            commonDivisor = std::gcd(commonDivisor, v);
        }
        if (commonDivisor == 0) {
            commonDivisor = 1;
        }

        bits = std::bit_width(largestValue / commonDivisor);
        for (auto v : args) {
            push_back(v);
        }
    }


    DenseVector(DenseVector const&) = default;
    DenseVector(DenseVector&&) noexcept = default;
    DenseVector& operator=(DenseVector const&) = default;
    DenseVector& operator=(DenseVector&&) noexcept = default;

    /** Reserve memory for a certain number of integers.
     *
     * \param s number of integer
     */
    void reserve(size_t s) {
        auto a = bits*s;

        data.reserve(a/64+((a % 64)>0?1:0));
    }


    /** Push a integer at the end of the vector
     *
     * \param value value to push onto the vector. Only the first `bits` are
     * being stored.
     */
    void push_back(uint64_t value) {
        assert(value % commonDivisor == 0);
        assert(value <= largestValue);
        value = value / commonDivisor;
        assert(std::bit_width(value) <= bits);
        auto empty = data.size()*64-bitCount;
        if (empty == 0) {
            data.push_back(value);
            bitCount += bits;
            return;
        }
        if (empty >= bits) {
            data.back() |= value << (64-empty);
            bitCount += bits;
            return;
        }

        data.back() |= value << (64-empty); // only pushes some bits, some are being dropped
        data.push_back(value >> empty);
        bitCount += bits;
    }

    /** Read integer at a certain position
     */
    auto operator[](size_t i) const -> uint64_t {
        return access(i);
    }

    /** Read integer at a certain position
     */
    auto access(size_t i) const -> uint64_t {
        auto begin = i * bits;
        auto end   = begin+bits-1; // The end is inclusive
        assert(begin == end || begin < bitCount);
        assert(end <= bitCount);
        assert(end >= begin);
        assert(end-begin <= 64);

        auto startI = begin/64;
        auto endI   = end/64;
        auto startOffset = begin % 64;

        auto mask = (uint64_t{1}<<(end-begin+size_t{1}))-1;

        auto value = [&]() -> uint64_t {
            // bits are all located inside a single uint64_t entry
            if (startI == endI) {
                return (data[startI] >> startOffset) & mask;
            }

            // bits are located in two adjecent uint64_t entries
            auto p1 = data[startI] >> startOffset;
            auto p2 = data[endI] << (64-startOffset);
            return (p1 | p2) & mask;

        }();

        return value * commonDivisor;
    }

    size_t size() const {
        return bitCount / bits; // Always a whole number
    }

    /** Number of bits required to store all integers
     */
    size_t bit_size() const {
        return bitCount;
    }

    /** Number of bits per entry
     */
    size_t entry_size() const {
        return bits;
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(data, bitCount, bits);
    }

};

}
