// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "bitset_popcount.h"
#include "concepts.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <libsais.h>
#include <libsais64.h>
#include <numeric>
#include <span>
#include <stdexcept>
#include <tuple>
#include <vector>

#if __has_include(<cereal/types/bitset.hpp>)
#include <cereal/types/bitset.hpp>
#endif


#if LIBSAIS_OPENMP
#   include <omp.h>
#endif

namespace fmindex_collection {

// helper to compute 'alignas' values
constexpr inline auto alignAsValue(size_t bits) -> size_t {
    if (bits % 512 == 0) return 64;
    if (bits % 256 == 0) return 32;
    if (bits % 128 == 0) return 16;
    if (bits %  64 == 0)  return 8;
    return 1;
}
// helper to compute 'alignas' values
constexpr inline auto minAlignAsValue(auto... bits) {
    return std::min(alignAsValue(bits)...);
}

template <size_t N, bool Align=true>
struct alignas(std::max(alignof(std::bitset<N>), Align?alignAsValue(N):size_t{1})) AlignedBitset {
    std::bitset<N> bits;

    decltype(auto) operator[](size_t i) {
        return bits[i];
    }
    decltype(auto) operator[](size_t i) const {
        return bits[i];
    }
    auto count() const {
        return bits.count();
    }

    template <typename Archive>
    void save(Archive& ar) const {
        saveBV(bits, ar);
    }

    template <typename Archive>
    void load(Archive& ar) {
        loadBV(bits, ar);
    }

};

inline auto createSA64(std::span<uint8_t const> input, size_t threadNbr) -> std::vector<uint64_t> {
    assert(uint64_t{input.size()} < std::numeric_limits<int64_t>::max());
    auto sa = std::vector<uint64_t>(input.size());
    if (input.size() == 0) {
        return sa;
    }
#if LIBSAIS_OPENMP
    auto r = libsais64_omp(input.data(), reinterpret_cast<int64_t*>(sa.data()), input.size(), 0, nullptr, threadNbr);
#else
    (void)threadNbr; // Unused if no openmp is available
    auto r = libsais64(input.data(), reinterpret_cast<int64_t*>(sa.data()), input.size(), 0, nullptr);
#endif

    if (r != 0) { throw std::runtime_error("something went wrong constructing the SA"); }
    return sa;
}

inline auto createSA32(std::span<uint8_t const> input, size_t threadNbr) -> std::vector<uint32_t> {
    assert(input.size() < std::numeric_limits<int32_t>::max());
    auto sa = std::vector<uint32_t>(input.size());
    if (input.size() == 0) {
        return sa;
    }
#if LIBSAIS_OPENMP
    auto r = libsais_omp(input.data(), reinterpret_cast<int32_t*>(sa.data()), input.size(), 0, nullptr, threadNbr);
#else
    (void)threadNbr; // Unused if no openmp is available
    auto r = libsais(input.data(), reinterpret_cast<int32_t*>(sa.data()), input.size(), 0, nullptr);
#endif

    if (r != 0) { throw std::runtime_error("something went wrong constructing the SA"); }
    return sa;
}

inline auto createBWT64(std::span<uint8_t const> input, std::span<uint64_t const> sa) -> std::vector<uint8_t> {
    assert(input.size() == sa.size());
    auto bwt = std::vector<uint8_t>{};
    bwt.resize(input.size());
    for (size_t i{0}; i < sa.size(); ++i) {
        bwt[i] = input[(sa[i] + input.size() - 1) % input.size()];
    }
    return bwt;
}

inline auto createBWT32(std::span<uint8_t const> input, std::span<uint32_t const> sa) -> std::vector<uint8_t> {
    assert(input.size() == sa.size());
    auto bwt = std::vector<uint8_t>{};
    bwt.resize(input.size());
    for (size_t i{0}; i < sa.size(); ++i) {
        bwt[i] = input[(sa[i] + input.size() - 1) % input.size()];
    }
    return bwt;
}

auto createSequences(Sequences auto const& _input, bool reverse=false) -> std::tuple<size_t, std::vector<uint8_t>, std::vector<size_t>> {
    // compute total numbers of bytes of the text including delimiters "$"
    size_t totalSize{};
    for (auto const& l : _input) {
        totalSize += l.size() + 1;
    }

    // our concatenated sequences with delimiters
    auto inputText = std::vector<uint8_t>{};
    inputText.reserve(totalSize);

    // list of sizes of the individual sequences
    auto inputSizes = std::vector<size_t>{};
    inputSizes.reserve(_input.size());

    for (auto const& l : _input) {
        if (not reverse) {
            inputText.insert(inputText.end(), begin(l), end(l));
        } else {
//!TODO hack for clang, broken in clang 15
#if __clang__
            auto l2 = std::vector<uint8_t>(l);
            std::ranges::reverse(l2);
#else
            auto l2 = std::views::reverse(l);
#endif
            inputText.insert(inputText.end(), begin(l2), end(l2));
        }

        // fill with delimiters/zeros
        inputText.resize(inputText.size() + 1);

        inputSizes.emplace_back(l.size()+1);
    }
    return {totalSize, inputText, inputSizes};
}

auto createSequencesAndReverse(Sequences auto const& _input) -> std::tuple<size_t, std::vector<uint8_t>, std::vector<size_t>> {
    // compute total numbers of bytes of the text including delimiters "$"
    size_t totalSize{};
    for (auto const& l : _input) {
        totalSize += l.size()+1;
    }
    totalSize = totalSize*2; // including reverse text

    // our concatenated sequences with delimiters
    auto inputText = std::vector<uint8_t>{};
    inputText.reserve(totalSize);

    // list of sizes of the individual sequences
    auto inputSizes = std::vector<size_t>{};
    inputSizes.reserve(_input.size());

    // add text
    for (auto const& l : _input) {
        inputText.insert(inputText.end(), begin(l), end(l));

        // fill with delimiters/zeros
        inputText.resize(inputText.size() + 1);

        inputSizes.emplace_back(l.size()+1);
    }

    // add reversed text
    for (auto const& l : std::views::reverse(_input)) {
        auto l2 = std::views::reverse(l);
        inputText.insert(inputText.end(), begin(l2), end(l2));

        // fill with delimiters/zeros
        inputText.resize(inputText.size() + 1);

        inputSizes.emplace_back(l.size()+1);
    }

    return {totalSize, inputText, inputSizes};
}



inline auto createSA_32(std::span<uint32_t const> input, size_t threadNbr) -> std::vector<int32_t> {
    //!TODO call to libsais_int_omp seems not correct
    auto sa = std::vector<int32_t>(input.size());
    if (input.size() == 0) {
        return sa;
    }
#if LIBSAIS_OPENMP
    auto r = libsais_int_omp((int32_t*)input.data(), sa.data(), input.size(), 65536, 0, threadNbr);
#else
    (void)threadNbr; // Unused if no openmp is available
    auto r = libsais_int((int32_t*)input.data(), sa.data(), input.size(), 65536, 0);
#endif

    if (r != 0) { throw std::runtime_error("something went wrong constructing the SA"); }
    return sa;
}


inline auto createBWT_32(std::span<uint32_t const> input, std::span<int32_t const> sa) -> std::vector<uint32_t> {
    assert(input.size() == sa.size());
    auto bwt = std::vector<uint32_t>{};
    bwt.resize(input.size());
    for (size_t i{0}; i < sa.size(); ++i) {
        bwt[i] = input[(sa[i] + input.size() - 1) % input.size()];
    }
    return bwt;
}

auto createSequences_32(Sequences auto const& _input, int samplingRate, bool reverse=false) -> std::tuple<size_t, std::vector<uint32_t>, std::vector<std::tuple<size_t, size_t>>> {
    // compute total numbers of bytes of the text including delimiters "$"
    size_t totalSize{};
    for (auto const& l : _input) {
        auto textLen    = l.size();
        auto delimLen = samplingRate - textLen % samplingRate; // Make sure it is always a multiple of samplingRate
        totalSize += textLen + delimLen;
    }

    // our concatenated sequences with delimiters
    auto inputText = std::vector<uint32_t>{};
    inputText.reserve(totalSize);

    // list of sizes of the individual sequences
    auto inputSizes = std::vector<std::tuple<size_t, size_t>>{};
    inputSizes.reserve(_input.size());


    for (auto const& l : _input) {
        auto ls = l.size();
        // number of delimiters ('$') which need to be added. It must be at least one, and it
        // has to make sure the text will be a multiple of samplingRate
        size_t delimCount = samplingRate - (ls % samplingRate);
        inputText.resize(inputText.size() + ls + delimCount, 0);

        if (not reverse) {
            std::ranges::copy(l, end(inputText) - ls - delimCount);
        } else {
//!TODO hack for clang, broken in clang 15
#if __clang__
            auto l2 = std::vector<uint32_t>(l);
            std::ranges::reverse(l2);
#else
            auto l2 = std::views::reverse(l);
#endif
            std::ranges::copy(l2, end(inputText) - ls - delimCount);
        }

        inputSizes.emplace_back(l.size(), delimCount);
    }
    return {totalSize, inputText, inputSizes};
}

struct IntIterator {
    size_t i;

    auto operator*() const -> size_t {
        return i;
    }

    auto operator++() -> IntIterator& {
        ++i;
        return *this;
    }

    auto operator<=>(IntIterator const&) const -> std::strong_ordering = default;
};

template <typename Index>
auto reconstructText(Index const& index, size_t seqNbr) -> std::vector<uint8_t> {
    auto r = std::vector<uint8_t>{};

    uint8_t c{};
    size_t idx = seqNbr;
    do {
        c = index.occ.symbol(idx);
        idx = index.occ.rank(idx, c);
        r.push_back(c);
    } while (c != 0);
    r.pop_back(); // remove last zero
    std::ranges::reverse(r);
    return r;
}

template <typename Index>
auto reconstructText(Index const& index) -> std::vector<std::vector<uint8_t>> {
    auto nbrOfSeq = index.occ.rank(index.size(), 0);
    auto texts = std::vector<std::vector<uint8_t>>{};
    auto seqIds = std::vector<std::tuple<size_t, size_t>>{};
    for (size_t i{}; i < nbrOfSeq; ++i) {
        texts.push_back(reconstructText(index, i));
        seqIds.emplace_back(std::get<0>(index.locate(i)), i);
    }
    std::ranges::sort(seqIds);
    auto res = std::vector<std::vector<uint8_t>>{};
    for (auto [seqId, idx] : seqIds) {
        res.emplace_back(std::move(texts[idx]));
    }
    return res;
}

namespace test {

/**
 * unittest helper function
 */
inline auto createInputData(std::vector<std::string> args) -> std::vector<std::vector<uint8_t>> {
    auto ret = std::vector<std::vector<uint8_t>>{};
    for (auto const& s : args) {
        ret.emplace_back();
        for (auto c : s) {
            ret.back().push_back(static_cast<uint8_t>(c));
        }
    }
    return ret;
}
}

}
