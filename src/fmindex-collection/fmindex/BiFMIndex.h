// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../string/concepts.h"
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

template <String_c String, SuffixArray_c TCSA = CSA>
struct BiFMIndex {
    static size_t constexpr Sigma = String::Sigma;

    String bwt;
    String bwtRev;
    std::array<size_t, Sigma+1> C{};
    TCSA   csa;

    BiFMIndex() = default;
    BiFMIndex(BiFMIndex&&) noexcept = default;

    BiFMIndex(std::span<uint8_t const> _bwt, std::span<uint8_t const> _bwtRev, TCSA _csa)
        : bwt{_bwt}
        , bwtRev{_bwtRev}
        , csa{std::move(_csa)}
    {
        for (auto c : _bwt) {
            C[c+1] += 1;
        }
        for (size_t i{1}; i < C.size(); ++i) {
            C[i] = C[i] + C[i-1];
        }

        assert(bwt.size() == bwtRev.size());
        if (bwt.size() != bwtRev.size()) {
            throw std::runtime_error("bwt don't have the same size: " + std::to_string(bwt.size()) + " " + std::to_string(bwtRev.size()));
        }
    }

    /**!\brief Creates a BiFMIndex with a specified sampling rate
     *
     * \param _input a list of sequences
     * \param samplingRate rate of the sampling
     */
    BiFMIndex(Sequences auto const& _input, size_t samplingRate, size_t threadNbr, bool useDelimiters = true, size_t seqOffset = 0) {

        auto [totalSize, inputText, inputSizes] = [&]() {
            if (useDelimiters) {
                return createSequences(_input);
            } else {
                return createSequencesWithoutDelimiter(_input);
            }
        }();

        //!TODO not properly pointing to the buffer, when reversing
        auto inputTextSpan = std::span<uint8_t const>{inputText};
        if (!useDelimiters) {
            // double text, so we get omage sorting
            auto halfSize = inputText.size();
            inputText.resize(halfSize*2);
            inputTextSpan = {inputText.begin(), inputText.begin()+halfSize};

            for (size_t i{0}; i < halfSize; ++i) {
                inputText[halfSize + i] = inputText[i];
            }
        }
        auto removeInvalidSA = [&](auto& sa) {
            if (!useDelimiters) { // using omega sorting, remove half of the entries
                auto halfSize = inputText.size() / 2;
                auto [first, last] = std::ranges::remove_if(sa, [&](auto e) {
                    return e >= halfSize;
                });
                sa.erase(first, last);
            }
        };

        if (totalSize < std::numeric_limits<int32_t>::max()) { // only 32bit SA required
            // create BurrowsWheelerTransform and CompressedSuffixArray
            auto [bwt, csa] = [&]() {
                auto sa  = createSA32(inputText, threadNbr);
                removeInvalidSA(sa);

                auto bwt = createBWT32(inputTextSpan, sa);
                auto csa = TCSA(std::move(sa), samplingRate, inputSizes, /*.reverse=*/false, /*.seqOffset=*/seqOffset);
                return std::make_tuple(std::move(bwt), std::move(csa));
            }();

            // create BurrowsWheelerTransform on reversed text
            auto bwtRev = [&]() {
                std::ranges::reverse(inputText);
                auto saRev  = createSA32(inputText, threadNbr);
                removeInvalidSA(saRev);

                auto bwtRev = createBWT32(inputTextSpan, saRev);
                return bwtRev;
            }();

            decltype(inputText){}.swap(inputText); // inputText memory can be deleted

            *this = BiFMIndex{bwt, bwtRev, std::move(csa)};
        } else { // required 64bit SA required
            // create BurrowsWheelerTransform and CompressedSuffixArray
            auto [bwt, csa] = [&]() {
                auto sa  = createSA64(inputText, threadNbr);
                removeInvalidSA(sa);
                auto bwt = createBWT64(inputTextSpan, sa);
                auto csa = TCSA(std::move(sa), samplingRate, inputSizes);
                return std::make_tuple(std::move(bwt), std::move(csa));
            }();

            // create BurrowsWheelerTransform on reversed text
            auto bwtRev = [&]() {
                std::ranges::reverse(inputText);
                auto saRev  = createSA64(inputText, threadNbr);
                removeInvalidSA(saRev);
                auto bwtRev = createBWT64(inputTextSpan, saRev);
                return bwtRev;
            }();

            decltype(inputText){}.swap(inputText); // inputText memory can be deleted

            *this = BiFMIndex{bwt, bwtRev, std::move(csa)};
        }
    }

    auto operator=(BiFMIndex const&) -> BiFMIndex& = delete;
    auto operator=(BiFMIndex&&) noexcept -> BiFMIndex& = default;

/*    size_t memoryUsage() const requires OccTableMemoryUsage<Table> {
        return occ.memoryUsage() + occRev.memoryUsage() + csa.memoryUsage();
    }*/

    size_t size() const {
        return bwt.size();
    }

    auto locate(size_t idx) const -> std::tuple<size_t, size_t> {
        if constexpr (requires(String t) {{ t.hasValue(size_t{}) }; }) {
            bool v = bwt.hasValue(idx);
            uint64_t steps{};
            while(!v) {
                idx = bwt.rank_symbol(idx);
                steps += 1;
                v = bwt.hasValue(idx);
            }
            auto [chr, pos] = csa.value(idx);
            return {chr, pos+steps};

        } else {
            auto opt = csa.value(idx);
            uint64_t steps{};
            while(!opt) {
                if constexpr (requires(String t) { { t.rank_symbol(size_t{}) }; }) {
                    idx = bwt.rank_symbol(idx);
                } else {
                    auto symb = bwt.symbol(idx);
                    idx = bwt.rank(idx, symb) + C[symb];
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
        ar(bwt, bwtRev, C, csa);
    }
};

}
#ifdef WORKAROUND_GCC_MACOS14
    #pragma GCC diagnostic pop
#endif
