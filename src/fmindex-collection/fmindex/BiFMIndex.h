// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../string/FlattenedBitvectors2L.h"
#include "../string/concepts.h"
#include "../suffixarray/SparseArray.h"
#include "../suffixarray/utils.h"
#include "../utils.h"

#include <algorithm>
#include <cassert>

namespace fmc {

template <size_t TSigma, template <size_t> typename String = string::FlattenedBitvectors_512_64k, SparseArray_c SparseArray = suffixarray::SparseArray<std::tuple<uint32_t, uint32_t>>, bool TDelim=true, bool TReuseRev=false>
    requires String_c<String<TSigma>>
struct BiFMIndex {
    using ADEntry = SparseArray::value_t;

    using NoDelim = BiFMIndex<TSigma, String, SparseArray, false, TReuseRev>;
    using ReuseRev = BiFMIndex<TSigma, String, SparseArray, TDelim, true>;

    static size_t constexpr Sigma     = TSigma;
    static size_t constexpr FirstSymb = TDelim?1:0;
    static bool constexpr Delim_v     = TDelim;
    static bool constexpr ReuseRev_v  = TReuseRev;

    // Set RevBwtType to std::nullptr_t to indicate that it should not be used
    using RevBwtType = std::conditional_t<TReuseRev, std::nullptr_t, String<Sigma>>;
    String<Sigma> bwt;
    RevBwtType bwtRev;
    std::array<size_t, Sigma+1> C{};
    SparseArray annotatedArray;

    BiFMIndex() = default;
    BiFMIndex(BiFMIndex&&) noexcept = default;

    BiFMIndex(std::span<uint8_t const> _bwt, std::span<uint8_t const> _bwtRev, SparseArray _annotatedArray)
        requires(!TReuseRev)
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

    BiFMIndex(std::span<uint8_t const> _bwt, SparseArray _annotatedArray)
        requires(TReuseRev)
        : bwt{_bwt}
        , C{computeC(bwt)}
        , annotatedArray{std::move(_annotatedArray)}
    {}


    /*
     * \param includeReversedInput assumes that the input data also has the reversed input data
     */
    BiFMIndex(Sequence auto const& _sequence, SparseArray const& _annotatedSequence, size_t _threadNbr, bool includeReversedInput = false) {
        if (_sequence.size() >= std::numeric_limits<size_t>::max()/2) {
            throw std::runtime_error{"sequence is longer than what this system is capable of handling"};
        }

        bool omegaSorting = !Delim_v; // Use omega sorting if no delimiter is being used

        // copy text into custom buffer
        auto inputText = createInputText(_sequence, omegaSorting, includeReversedInput);

        // create bwt, bwtRev and annotatedArray
        auto [_bwt, _annotatedArray] = createBWTAndAnnotatedArray(inputText, _annotatedSequence, _threadNbr, omegaSorting);


        if constexpr (!TReuseRev) {
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
            bwtRev = {_bwtRev};
        }

        // initialize this BiFMIndex properly
        bwt = {_bwt};
        C = computeC(bwt);
        annotatedArray = std::move(_annotatedArray);
    }


    /**!\brief Creates a BiFMIndex with a specified sampling rate
     *
     * \param _input a list of sequences
     * \param samplingRate rate of the sampling
     * \param includeReversedInput also adds all input and their reversed text
     */
    BiFMIndex(Sequences auto const& _input, size_t samplingRate, size_t threadNbr, size_t seqOffset = 0, bool includeReversedInput = false) {
        auto [totalSize, inputText, inputSizes] = createSequences(_input, /*._addReversed=*/includeReversedInput, /*._useDelimiters=*/Delim_v);

        size_t refId{0};
        size_t pos{0};

        //!TODO what about empty strings?
        while (inputSizes[refId] == pos) {
            refId += 1;
            assert(refId < inputSizes.size());
        }
        auto const startRefId = refId;


        auto annotatedSequence = SparseArray {
            std::views::iota(size_t{0}, totalSize) | std::views::transform([&](size_t phase) -> std::optional<ADEntry> {
                if (phase == 0) { // restarting
                    refId = startRefId;
                    pos = 0;
                }
                if (!includeReversedInput || phase*2 < totalSize) { // going forward
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
                } else { // going backwards
                    assert(refId < inputSizes.size());
                    assert(pos < inputSizes[refId]);

                    auto ret = std::optional<ADEntry>{std::nullopt};

                    if (pos % samplingRate == 0) {
                        auto _refId = _input.size() + inputSizes.size() - refId-1+seqOffset;
                        size_t extra = Delim_v?1:0;
                        auto _pos   = (inputSizes[refId] - pos + inputSizes[refId] - 1 - extra) % inputSizes[refId];
                        ret = std::make_tuple(_refId, _pos);
                    }

                    ++pos;
                    if (inputSizes[refId] == pos) {
                        refId += 1;
                        pos = 0;
                    }
                    return ret;
                }
            })
        };

        *this = BiFMIndex{inputText, annotatedSequence, threadNbr, /*includeReversedInput=*/false};
    }

    auto operator=(BiFMIndex const&) -> BiFMIndex& = delete;
    auto operator=(BiFMIndex&& _other) noexcept -> BiFMIndex& = default;

    size_t size() const {
        return bwt.size();
    }

    using LEntry = decltype(std::tuple_cat(ADEntry{}, std::tuple<size_t>{}));
    auto locate(size_t idx) const -> LEntry {
        if constexpr (requires(String<Sigma> t) {{ t.hasValue(size_t{}) }; }) {
            bool v = bwt.hasValue(idx);
            uint64_t steps{};
            while (!v) {
                idx = bwt.rank_symbol(idx);
                steps += 1;
                v = bwt.hasValue(idx);
            }
            return std::tuple_cat(*annotatedArray.value(idx), std::tuple<size_t>{steps});
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
            return std::tuple_cat(*opt, std::tuple<size_t>{steps});
        }
    }

    auto single_locate_step(size_t idx) const -> std::optional<ADEntry> {
        return annotatedArray.value(idx);
    }


    template <typename Archive>
    void serialize(this auto&& self, Archive& ar) {
        ar(self.bwt, self.C, self.annotatedArray);
        if constexpr (!std::same_as<RevBwtType, std::nullptr_t>) {
            ar(self.bwtRev);
        }
    }
};

}
