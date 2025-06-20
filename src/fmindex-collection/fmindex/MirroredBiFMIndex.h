// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../string/concepts.h"
#include "../suffixarray/CSA.h"
#include "../utils.h"

#include <algorithm>

namespace fmc {

template <String_c String, SuffixArray_c TCSA = CSA>
struct MirroredBiFMIndex {
    using ADEntry = std::tuple<size_t, size_t>;
    static size_t constexpr Sigma = String::Sigma;

    String bwt;
    std::array<size_t, Sigma+1> C{};
    TCSA   csa;

    MirroredBiFMIndex() = default;
    MirroredBiFMIndex(std::span<uint8_t const> _bwt, TCSA _csa)
        : bwt{_bwt}
        , csa{std::move(_csa)}
    {
        for (auto c : _bwt) {
            C[c+1] += 1;
        }
        for (size_t i{1}; i < C.size(); ++i) {
            C[i] = C[i] + C[i-1];
        }
    }

    /**!\brief Creates a MirroredBiFMIndex with a specified sampling rate
     *
     * \param _input a list of sequences
     * \param samplingRate rate of the sampling
     */
    MirroredBiFMIndex(Sequences auto const& _input, size_t samplingRate, size_t threadNbr, bool useDelimiters = true) {
        auto [totalSize, inputText, inputSizes] = [&]() {
            if (useDelimiters) {
                return createSequencesAndReverse(_input);
            } else {
                return createSequencesAndReverseWithoutDelimiter(_input);
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


        // Check only valid characters are used
        assert([&]() {
            for (auto c : inputText) {
                if (c >= Sigma) return false;
            }
            return true;
        }());

        // create BurrowsWheelerTransform and CompressedSuffixArray
        auto [bwt, csa] = [&]() {
            auto sa  = createSA64(inputText, threadNbr);
            removeInvalidSA(sa);

            auto bwt = createBWT64(inputTextSpan, sa);
            auto csa = TCSA(std::move(sa), samplingRate, inputSizes);
            return std::make_tuple(std::move(bwt), std::move(csa));
        }();

        decltype(inputText){}.swap(inputText); // inputText memory can be deleted

        *this = MirroredBiFMIndex{bwt, std::move(csa)};
    }

    size_t size() const {
        return bwt.size();
    }

    auto locate(size_t idx) const -> std::tuple<ADEntry, size_t> {
        if constexpr (requires(String t) {{ t.hasValue(size_t{}) }; }) {
            bool v = bwt.hasValue(idx);
            uint64_t steps{};
            while(!v) {
                idx = bwt.rank_symbol(idx);
                steps += 1;
                v = bwt.hasValue(idx);
            }
            return {csa.value(idx), steps};
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
            return {*opt, steps};
        }
    }

    auto single_locate_step(size_t idx) const -> std::optional<ADEntry> {
        return csa.value(idx);
    }


    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bwt, C, csa);
    }
};

}
