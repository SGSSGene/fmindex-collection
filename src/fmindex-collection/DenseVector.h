// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "cereal_tag.h"

#include <cassert>
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
struct DenseVector {
    std::vector<uint64_t> data; // buffer where the data is being stored
    size_t bitCount{};          // numbers of used bits
    size_t bits{};              // number of bits per entry

    /**
     * C'Tor for serialization purposes
     */
    DenseVector(cereal_tag) {}

    /**
     * \param _bits number of bits used for each integer entry
     */
    DenseVector(size_t _bits)
        : bits{_bits}
    {}

    /**
     * Create a DenseVector and auto determine the number of bits required
     */
    DenseVector(std::initializer_list<uint64_t> args) {
        for (auto v : args) {
            size_t ct = 0;
            while (v > 0) {
                v = v >> 1;
                ct += 1;
            }
            bits = std::max(bits, ct);
        }
        for (auto v : args) {
            push_back(v);
        }
    }

    DenseVector(DenseVector&&) noexcept = default;
    DenseVector(DenseVector const&) = default;
    DenseVector& operator=(DenseVector&&) noexcept = default;
    DenseVector& operator=(DenseVector const&) = default;

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
    auto operator[](size_t i) const {
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
//        auto endOffset   = end % 64;

        auto mask = (1ull<<(end-begin+1ull))-1ull;
        if (startI == endI) {
            return (data[startI] >> startOffset) & mask;
        }

        auto p1 = data[startI] >> startOffset;
        auto p2 = data[endI] << (64-startOffset);
        return (p1 | p2) & mask;
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
