// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../occtable/concepts.h"
#include "../suffixarray/CSA.h"
#include "../utils.h"

#include <algorithm>

#define GCC_COMPILER (defined(__GNUC__) && !defined(__clang__))

#if defined(__GNUC__) && !defined(__clang__) && defined(__MACH__)
    #define WORKAROUND_GCC_MACOS14
#else
#endif

//!WORKAROUND: only triggers on macos-14 with gcc
#ifdef WORKAROUND_GCC_MACOS14
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif



namespace fmindex_collection {

template <OccTable Table, SuffixArray_c TCSA = CSA>
struct BiFMIndex {
    static size_t constexpr Sigma = Table::Sigma;

    Table  occ;
    Table  occRev;
    TCSA   csa;

    BiFMIndex() = default;
    BiFMIndex(BiFMIndex&&) noexcept = default;

    BiFMIndex(std::span<uint8_t const> bwt, std::span<uint8_t const> bwtRev, TCSA _csa)
        : occ{bwt}
        , occRev{bwtRev}
        , csa{std::move(_csa)}
    {
        assert(bwt.size() == bwtRev.size());
        assert(occ.size() == occRev.size());
        if (bwt.size() != bwtRev.size()) {
            throw std::runtime_error("bwt don't have the same size: " + std::to_string(bwt.size()) + " " + std::to_string(bwtRev.size()));
        }
        if (occ.size() != occRev.size()) {
            throw std::runtime_error("occ don't have the same size: " + std::to_string(occ.size()) + " " + std::to_string(occRev.size()));
        }
        // compute last row
        auto ct = std::array<uint64_t, Sigma>{};
        for (auto v : bwt) {
            ct[v] += 1;
        }
        for (size_t i{1}; i < ct.size(); ++i) {
            ct[i] = ct[i-1] + ct[i];
        }
        // check last row is correct
        for (size_t sym{0}; sym < Sigma; ++sym) {
            if (occ.rank(occ.size(), sym) != ct[sym]) {
                auto e = std::string{"Wrong rank for the last entry."}
                    + " Got different values for forward index."
                    + " sym: " + std::to_string(sym)
                    + " got: " + std::to_string(occ.rank(occ.size(), sym))
                    + " expected: " + std::to_string(ct[sym]);
                throw std::runtime_error(e);
            }
            if (occRev.rank(occRev.size(), sym) != ct[sym]) {
                auto e = std::string{"Wrong rank for the last entry."}
                    + " Got different values for reverse index."
                    + " sym: " + std::to_string(sym)
                    + " got: " + std::to_string(occRev.rank(occRev.size(), sym))
                    + " expected: " + std::to_string(ct[sym]);
                throw std::runtime_error(e);
            }
        }
        if constexpr (requires(Table t) {{ t.hasValue(size_t{}) }; }) {
            for (size_t i{0}; i < occ.size(); ++i) {
                if (csa.value(i).has_value()) {
                    occ.setValue(i);
                }
            }
        }
    }

    /**!\brief Creates a BiFMIndex with a specified sampling rate
     *
     * \param _input a list of sequences
     * \param samplingRate rate of the sampling
     */
    BiFMIndex(Sequences auto const& _input, size_t samplingRate, size_t threadNbr) {
        auto [totalSize, inputText, inputSizes] = createSequences(_input);

        if (totalSize < std::numeric_limits<int32_t>::max()) { // only 32bit SA required
            // create BurrowsWheelerTransform and CompressedSuffixArray
            auto [bwt, csa] = [&]() {
                auto sa  = createSA32(inputText, threadNbr);
                auto bwt = createBWT32(inputText, sa);
                auto csa = TCSA(std::move(sa), samplingRate, inputSizes);
                return std::make_tuple(std::move(bwt), std::move(csa));
            }();

            // create BurrowsWheelerTransform on reversed text
            auto bwtRev = [&]() {
                std::ranges::reverse(inputText);
                auto saRev  = createSA32(inputText, threadNbr);
                auto bwtRev = createBWT32(inputText, saRev);
                return bwtRev;
            }();

            decltype(inputText){}.swap(inputText); // inputText memory can be deleted

            *this = BiFMIndex{bwt, bwtRev, std::move(csa)};
        } else { // required 64bit SA required
            // create BurrowsWheelerTransform and CompressedSuffixArray
            auto [bwt, csa] = [&]() {
                auto sa  = createSA64(inputText, threadNbr);
                auto bwt = createBWT64(inputText, sa);
                auto csa = TCSA(std::move(sa), samplingRate, inputSizes);
                return std::make_tuple(std::move(bwt), std::move(csa));
            }();

            // create BurrowsWheelerTransform on reversed text
            auto bwtRev = [&]() {
                std::ranges::reverse(inputText);
                auto saRev  = createSA64(inputText, threadNbr);
                auto bwtRev = createBWT64(inputText, saRev);
                return bwtRev;
            }();

            decltype(inputText){}.swap(inputText); // inputText memory can be deleted

            *this = BiFMIndex{bwt, bwtRev, std::move(csa)};
        }
    }

    auto operator=(BiFMIndex const&) -> BiFMIndex& = delete;
    auto operator=(BiFMIndex&&) noexcept -> BiFMIndex& = default;

    size_t memoryUsage() const requires OccTableMemoryUsage<Table> {
        return occ.memoryUsage() + occRev.memoryUsage() + csa.memoryUsage();
    }

    size_t size() const {
        return occ.size();
    }

    auto locate(size_t idx) const -> std::tuple<size_t, size_t> {
        if constexpr (requires(Table t) {{ t.hasValue(size_t{}) }; }) {
            bool v = occ.hasValue(idx);
            uint64_t steps{};
            while(!v) {
                idx = occ.rank_symbol(idx);
                steps += 1;
                v = occ.hasValue(idx);
            }
            auto [chr, pos] = csa.value(idx);
            return {chr, pos+steps};

        } else {
            auto opt = csa.value(idx);
            uint64_t steps{};
            while(!opt) {
                if constexpr (requires(Table t) { { t.rank_symbol(size_t{}) }; }) {
                    idx = occ.rank_symbol(idx);
                } else {
                    idx = occ.rank(idx, occ.symbol(idx));
                }
                steps += 1;
                opt = csa.value(idx);
            }
            auto [chr, pos] = *opt;
            return {chr, pos+steps};
        }
    }

    auto single_locate_step(size_t idx) const -> std::optional<std::tuple<size_t, size_t>> {
        return csa.value(idx);
    }


    template <typename Archive>
    void serialize(Archive& ar) {
        ar(occ, occRev, csa);
    }
};

}
#ifdef WORKAROUND_GCC_MACOS14
    #pragma GCC diagnostic pop
#endif
