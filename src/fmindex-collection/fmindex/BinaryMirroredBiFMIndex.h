// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../bitvector/concepts.h"
#include "../suffixarray/CSA.h"
#include "../utils.h"

#include <algorithm>

namespace fmc {

template <Bitvector_c Vector, SuffixArray_c TCSA = CSA>
struct BinaryMirroredBiFMIndex {
    using ADEntry = std::tuple<size_t, size_t>;
    static size_t constexpr Sigma = 2;

    Vector bwt;
    std::array<size_t, Sigma+1> C{};
    TCSA   csa;

    BinaryMirroredBiFMIndex() = default;
    BinaryMirroredBiFMIndex(std::span<uint8_t const> _bwt, TCSA _csa)
        : bwt{_bwt}
        , csa{std::move(_csa)}
    {
        for (auto c : _bwt) {
            assert(c < Sigma); // must fit into the alphabet
            C[c+1] += 1;
        }
        for (size_t i{1}; i < C.size(); ++i) {
            C[i] = C[i] + C[i-1];
        }
    }

    /**!\brief Creates a BinaryMirroredBiFMIndex with a specified sampling rate
     *
     * \param _input a list of sequences
     * \param samplingRate rate of the sampling
     */
    BinaryMirroredBiFMIndex(Sequences auto const& _input, size_t samplingRate, size_t threadNbr, bool useDelimiters = true) {
        auto [totalSize, inputText, inputSizes] = [&]() {
            if (useDelimiters) {
                return createSequencesAndReverse(_input);
            } else {
                return createSequencesAndReverseWithoutDelimiter(_input);
            }
        }();

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
            auto bwt = createBWT64(inputText, sa);
            auto csa = TCSA(std::move(sa), samplingRate, inputSizes);
            return std::make_tuple(std::move(bwt), std::move(csa));
        }();

        decltype(inputText){}.swap(inputText); // inputText memory can be deleted

        *this = BinaryMirroredBiFMIndex{bwt, std::move(csa)};
    }

    size_t size() const {
        return bwt.size();
    }

    auto locate(size_t idx) const -> std::tuple<ADEntry, size_t> {
        if constexpr (requires(Vector t) {{ t.hasValue(size_t{}) }; }) {
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
                if constexpr (requires(Vector t) { { t.rank_symbol(size_t{}) }; }) {
                    idx = bwt.rank_symbol(idx);
                } else {
                    auto symb = bwt.symbol(idx);
                    if (symb) {
                        idx = bwt.rank(idx) + C[symb] ;
                    } else {
                        idx = idx - bwt.rank(idx) + C[symb];
                    }
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
