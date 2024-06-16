// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../rankvector/concepts.h"
#include "concepts.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <span>

#if __has_include(<cereal/archives/binary.hpp>)
#include <cereal/types/array.hpp>
#endif

/**
 * A run length encoded bwt, see unpublished paper, but with recursive bitvectors
 *
 */
namespace fmindex_collection::occtable {

/**
 * Literal class type that wraps a constant expression string.
 *
 * Uses implicit conversion to allow templates to *seemingly* accept constant strings.
 */
template<size_t N>
struct StringLiteral {
    constexpr StringLiteral(char const (&str)[N]) {
        std::copy_n(str, N, value.data());
    }
    std::array<char, N> value;

    operator std::string() const {
        return std::string{value.data()};
    }
};

template <RankVector Vector, StringLiteral Name, StringLiteral Extension>
struct GenericOccTable {
    static constexpr size_t Sigma = Vector::Sigma;

    Vector                    vector;
    std::array<size_t, Sigma> C{};

    GenericOccTable() = default;
    GenericOccTable(std::span<uint8_t const> _symbols)
        : vector{_symbols}
    {
        for (auto c : _symbols) {
            C[c+1] += 1;
        }
        for (size_t i{1}; i < C.size(); ++i) {
            C[i] = C[i] + C[i-1];
        }
    }

    auto prefetch(uint64_t idx) const {
        (void)idx;
        vector.prefetch();
    }

    size_t size() const {
        return vector.size();
    }

    uint8_t symbol(size_t idx) const {
        return vector.symbol(idx);
    }

    uint64_t rank(size_t idx, uint64_t symb) const {
        return vector.rank(idx, symb) + C[symb];
    }

    uint64_t prefix_rank(size_t idx, uint64_t symb) const {
        return vector.prefix_rank(idx, symb);
    }

    auto all_ranks(size_t idx) const -> std::tuple<std::array<uint64_t, Sigma>, std::array<uint64_t, Sigma>> {
        auto [rs, prs] = vector.all_ranks_and_prefix_ranks(idx);
        for (size_t i{0}; i < rs.size(); ++i) {
            rs[i] += C[i];
        }
        return {rs, prs};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(vector, C);
    }

    static auto name() -> std::string {
        return Name;
    }

    static auto extension() -> std::string {
        return Extension;
    }

};

}
