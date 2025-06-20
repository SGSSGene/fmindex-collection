// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../string/FlattenedBitvectors_L0L1.h"
#include "../string/concepts.h"
#include "../suffixarray/SparseArray.h"
#include "../suffixarray/utils.h"
#include "../utils.h"

#include <algorithm>
#include <cassert>

namespace fmc {

template <size_t TSigma, template <size_t> typename String = string::FlattenedBitvectors_512_64k, SparseArray_c SparseArray = suffixarray::SparseArray<std::tuple<size_t, size_t>>>
    requires String_c<String<TSigma>>
struct BiFMIndex {
    using ADEntry = SparseArray::value_t;

    static size_t constexpr Sigma     = TSigma;
    static size_t constexpr FirstSymb = 1;

    String<Sigma> bwt;
    String<Sigma> bwtRev;
    std::array<size_t, Sigma+1> C{};
    SparseArray annotatedArray;

    BiFMIndex() = default;
    BiFMIndex(BiFMIndex&&) noexcept = default;

    BiFMIndex(std::span<uint8_t const> _bwt, std::span<uint8_t const> _bwtRev, SparseArray _annotatedArray)
        : bwt{_bwt}
        , bwtRev{_bwtRev}
        , C{computeC(bwt)}
        , annotatedArray{std::move(_annotatedArray)}
    {
        assert(bwt.size() == bwtRev.size());
        if (bwt.size() != bwtRev.size()) {
            throw std::runtime_error("bwt don't have the same size: " + std::to_string(bwt.size()) + " " + std::to_string(bwtRev.size()));
        }
    }

    BiFMIndex(Sequence auto const& _sequence, SparseArray const& _annotatedSequence, size_t _threadNbr, bool omegaSorting = false) {
        if (_sequence.size() >= std::numeric_limits<size_t>::max()/2) {
            throw std::runtime_error{"sequence is longer than what this system is capable of handling"};
        }

        // copy text into custom buffer
        auto inputText = createInputText(_sequence, omegaSorting);

        // create bwt, bwtRev and annotatedArray
        auto [_bwt, _annotatedArray] = createBWTAndAnnotatedArray(inputText, _annotatedSequence, _threadNbr, omegaSorting);
        #if defined(__GNUC__) && !defined(__clang__)
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wstringop-overflow"
        std::ranges::reverse(inputText);
        #pragma GCC diagnostic pop

        #else
        std::ranges::reverse(inputText);

        #endif
        auto _bwtRev = createBWT(inputText, _threadNbr, omegaSorting);
        decltype(inputText){}.swap(inputText); // inputText memory can be deleted

        // initialize this BiFMIndex properly
        bwt = {_bwt};
        bwtRev = {_bwtRev};
        C = computeC(bwt);
        annotatedArray = std::move(_annotatedArray);
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

        size_t refId{0};
        size_t pos{0};

        //!TODO what about empty strings?
        while (inputSizes[refId] == pos) {
            refId += 1;
            assert(inputSizes.size() < refId);
        }

        auto annotatedSequence = SparseArray {
            std::views::iota(size_t{0}, totalSize) | std::views::transform([&](size_t) -> std::optional<ADEntry> {
                assert(refId < inputSizes.size());
                assert(pos < inputSizes[refId]);

                auto ret = std::optional<ADEntry>{std::nullopt};

                if (pos % samplingRate == 0) {
                    ret = std::make_tuple(refId+seqOffset, pos);
                }

                ++pos;
                if (inputSizes[refId] == pos) {
                    refId += 1;
                    pos = 0;
                }
                return ret;
            })
        };

        *this = BiFMIndex{inputText, annotatedSequence, threadNbr, !useDelimiters};
    }

    auto operator=(BiFMIndex const&) -> BiFMIndex& = delete;
    auto operator=(BiFMIndex&&) noexcept -> BiFMIndex& = default;

    size_t size() const {
        return bwt.size();
    }

    auto locate(size_t idx) const -> std::tuple<ADEntry, size_t> {
        if constexpr (requires(String<Sigma> t) {{ t.hasValue(size_t{}) }; }) {
            bool v = bwt.hasValue(idx);
            uint64_t steps{};
            while (!v) {
                idx = bwt.rank_symbol(idx);
                steps += 1;
                v = bwt.hasValue(idx);
            }
            return {*annotatedArray.value(idx), steps};

        } else {
            auto opt = annotatedArray.value(idx);
            uint64_t steps{};
            while (!opt) {
                if constexpr (requires(String<Sigma> t) { { t.rank_symbol(size_t{}) }; }) {
                    idx = bwt.rank_symbol(idx);
                } else {
                    auto symb = bwt.symbol(idx);
                    idx = bwt.rank(idx, symb) + C[symb];
                }
                steps += 1;
                opt = annotatedArray.value(idx);
            }
            return {*opt, steps};
        }
    }

    auto single_locate_step(size_t idx) const -> std::optional<ADEntry> {
        return annotatedArray.value(idx);
    }


    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bwt, bwtRev, C, annotatedArray);
    }
};

}
