// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <tuple>
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
template <size_t WordWidth, size_t MaxValue = std::numeric_limits<size_t>::max(), size_t EntriesPerMWord = 64/WordWidth>
struct FixedSuccinctVector2 {
    std::vector<uint64_t> data; // buffer where the data is being stored
    size_t entryCount{};        // number of entries
    constexpr static uint64_t mask{(~uint64_t{0}) >> (64-WordWidth)};

    FixedSuccinctVector2() = default;
    FixedSuccinctVector2(std::initializer_list<uint64_t> list) {
        for (auto s : list) {
            push_back(s);
        }
    }
    FixedSuccinctVector2(FixedSuccinctVector2 const&) = default;
    FixedSuccinctVector2(FixedSuccinctVector2&&) noexcept = default;
    FixedSuccinctVector2& operator=(FixedSuccinctVector2 const&) = default;
    FixedSuccinctVector2& operator=(FixedSuccinctVector2&&) noexcept = default;

    /** Reserve memory for a certain number of integers.
     *
     * \param s number of integer
     */
    void reserve(size_t s) {
        auto a = (s+EntriesPerMWord-1)/EntriesPerMWord;

        data.reserve(a);
    }

    void resize(size_t s) {
        entryCount = s;
        auto a = (s+EntriesPerMWord-1)/EntriesPerMWord;
        data.resize(a);
    }


    /** Push a integer at the end of the vector
     *
     * \param value value to push onto the vector. Only the first `WordWidth` are
     * being stored.
     */
    void push_back(uint64_t value) {
        assert(value <= MaxValue);
        assert(std::log2(value) < WordWidth);
        auto count = data.size()*EntriesPerMWord - entryCount;
        entryCount += 1;
        if (count == 0) {
            data.push_back(value);
            count = EntriesPerMWord;
        } else {
            set(entryCount-1, value);
        }
        assert(back() == value);
    }

    void emplace_back(uint64_t value) {
        push_back(value);
    }

    void set(size_t idx, uint64_t value) {
        assert(value <= MaxValue);
        assert(std::log2(value) < WordWidth);

        size_t aid    = idx / EntriesPerMWord;
        size_t lid    = idx % EntriesPerMWord;
        size_t offset = lid * WordWidth;

        data[aid] = (data[aid] &~ (mask << offset)) | (value << offset);
        assert(at(idx) == value);
    }

    /** Read integer at a certain position
     */
    auto operator[](size_t i) const {
        return at(i);
    }

    /** Read integer at a certain position
     */
    auto at(size_t idx) const -> uint64_t {
        size_t aid    = idx / EntriesPerMWord;
        size_t lid    = idx % EntriesPerMWord;
        size_t offset = lid*WordWidth;

        return (data[aid] >> offset) & mask;
    }

    void pop_back() {
        assert(entryCount > 0);
        entryCount -= 1;
        auto count = data.size()*EntriesPerMWord - entryCount;
        assert(count <= EntriesPerMWord);
        if (count == EntriesPerMWord) {
            data.pop_back();
        }
    }

    uint64_t back() const {
        return at(size()-1);
    }

    size_t size() const {
        return entryCount;
    }

    /** Number of bits required to store all integers
     */
    size_t bit_size() const {
        return data.size() * 64;
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(data, entryCount);
    }
};

}
