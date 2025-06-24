// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../string/concepts.h"
#include "../string/utils.h"
#include "../string/FlattenedBitvectors_L0L1.h"
#include "../suffixarray/SparseArray.h"
#include "../utils.h"

namespace fmc {

template <size_t TSigma, template <size_t> typename String = string::FlattenedBitvectors_512_64k, SparseArray_c SparseArray = suffixarray::SparseArray<std::tuple<uint32_t, uint32_t>>>
    requires String_c<String<TSigma>>
struct FMIndex {
    using ADEntry = SparseArray::value_t;

    static size_t constexpr Sigma = TSigma;

    String<Sigma>               bwt;
    std::array<size_t, Sigma+1> C{0};
    SparseArray annotatedArray;

    FMIndex() = default;
    FMIndex(FMIndex const&) = delete;
    FMIndex(FMIndex&&) noexcept = default;
    FMIndex(std::span<uint8_t const> _bwt, SparseArray _annotatedArray)
        : bwt{_bwt}
        , C{computeC(bwt)}
        , annotatedArray{std::move(_annotatedArray)}
    {}

    FMIndex(Sequence auto const& _sequence, SparseArray const& _annotatedSequence, size_t _threadNbr, bool omegaSorting = false) {
        // copy text into custom buffer
        auto inputText = createInputText(_sequence, omegaSorting);

        // create bwt, bwtRev and annotatedArray
        auto [_bwt, _annotatedArray] = createBWTAndAnnotatedArray(inputText, _annotatedSequence, _threadNbr, omegaSorting);
        decltype(inputText){}.swap(inputText); // inputText memory can be deleted

        // initialize this FMIndex properly
        bwt = {_bwt};
        C = computeC(bwt);
        annotatedArray = std::move(_annotatedArray);
    }

    FMIndex(Sequence auto const& _input, size_t samplingRate, size_t threadNbr, bool useDelimiters = true, size_t seqOffset = 0)
        : FMIndex{std::vector<std::vector<uint8_t>>{_input}, samplingRate, threadNbr, useDelimiters, seqOffset}
    {}


    /**!\brief Creates a FMIndex with a specified sampling rate
     *
     * \param _input a list of sequences
     * \param samplingRate rate of the sampling
     */
    FMIndex(Sequences auto const& _input, size_t samplingRate, size_t threadNbr, bool useDelimiters = true, size_t seqOffset = 0) {

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

        *this = FMIndex{inputText, annotatedSequence, threadNbr, !useDelimiters};
    }
    auto operator=(FMIndex const&) -> FMIndex& = delete;
    auto operator=(FMIndex&&) noexcept -> FMIndex& = default;


    size_t size() const {
        return bwt.size();
    }

    auto locate(size_t idx) const -> std::tuple<ADEntry, size_t> {
        auto opt = annotatedArray.value(idx);
        size_t steps{};
        while(!opt) {
            auto symb = bwt.symbol(idx);
            idx = bwt.rank(idx, symb) + C[symb];
            steps += 1;
            opt = annotatedArray.value(idx);
        }
        return {*opt, steps};
    }

    auto single_locate_step(size_t idx) const -> std::optional<ADEntry> {
        return annotatedArray.value(idx);
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bwt, C, annotatedArray);
    }
};

}
