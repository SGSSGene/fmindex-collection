// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../string/PairedFlattenedBitvectors2L_b.h"
#include "../string/AdapterStringKStep.h"
#include "../string/concepts.h"
#include "../suffixarray/SparseArray.h"
#include "../suffixarray/utils.h"
#include "../VectorBool.h"
#include "../utils.h"

#include <algorithm>
#include <cassert>

namespace fmc {

template <size_t Base>
constexpr static size_t my_pow_base(size_t e) {
    if (e == 0) return 1;
    if (e == 1) return Base;
    auto v = my_pow_base<Base>(e/2);
    v = v*v;
    if (e%2 == 1) {
        v *= Base;
    }
    return v;
}


template <size_t TSigma, template <size_t> typename TString = string::PairedFlattenedBitvectors_b_512_64k, SparseArray_c TSparseArray = suffixarray::SparseArray<std::tuple<uint32_t, uint32_t>>, bool TDelim=true, bool TReuseRev=false, size_t TKStep=2>
    requires String_c<TString<TSigma>>
struct BiFMIndexKStep {
    using DocumentEntries = TSparseArray::value_t;
    using ADEntry = DocumentEntries;
    using SparseArray = TSparseArray;

    using NoDelim = BiFMIndexKStep<TSigma, TString, TSparseArray, false, TReuseRev, TKStep>;
    using ReuseRev = BiFMIndexKStep<TSigma, TString, TSparseArray, TDelim, true, TKStep>;

    template <size_t TNewKStep>
    using SetKStep = BiFMIndexKStep<TSigma, TString, TSparseArray, TDelim, TReuseRev, TNewKStep>;

    static size_t constexpr Sigma     = TSigma;
    static size_t constexpr FirstSymb = TDelim?1:0;
    static bool constexpr Delim_v     = TDelim;
    static bool constexpr ReuseRev_v  = TReuseRev;

    static size_t constexpr KStep = TKStep;

    static constexpr auto bitct = std::bit_width(TSigma-1);
    static size_t constexpr SigmaKStep = size_t{1}<<(bitct*KStep);

    // If String class does not fulfill StringKStep_c concept, wrap it in
    // AdapterStringKStep such that it does.
    using StringKStep = std::conditional_t<
        StringKStep_c<TString<SigmaKStep>>,
        TString<SigmaKStep>,
        fmc::string::AdapterStringKStep<SigmaKStep, bitct, TString>
    >;

    using RevBwtKStepType = std::conditional_t<TReuseRev, std::nullptr_t, StringKStep>;

    StringKStep bwt_kstep;
    RevBwtKStepType bwtRev_kstep;

    std::array<size_t, Sigma+1> C{};
    std::array<size_t, SigmaKStep+1> C_kstep{};
    std::array<size_t, SigmaKStep+1> CRev_kstep{};
    SparseArray annotatedArray;
    VectorBool annotatedArrayIsKStep;

    BiFMIndexKStep() = default;
    BiFMIndexKStep(BiFMIndexKStep&&) noexcept = default;


    static auto helperSwapC(std::array<size_t, SigmaKStep+1> a) {
        auto values = std::array<size_t, KStep>{};
        auto increment = std::function<void(size_t pos)>{};

        increment = [&](size_t pos) {
            if (pos == KStep) return;
            values[pos] += 1;
            if (values[pos] == Sigma) {
                values[pos] = 0;
                increment(pos+1);
            }
        };
        auto asValueFwd = [](std::array<size_t, KStep> const& v) {
            size_t value{};
            for (size_t i{0}; i < KStep; ++i) {
                value = value * Sigma + v[i];
            }
            return value;
        };
        auto asValueBwd = [](std::array<size_t, KStep> const& v) {
            size_t value{};
            for (size_t i{0}; i < KStep; ++i) {
                value = value * Sigma + v[KStep -1 -i];
            }
            return value;
        };

        for (size_t i{0}; i < std::pow(Sigma, KStep); ++i) {
            auto pos1 = asValueFwd(values);
            auto pos2 = asValueBwd(values);
            if (pos1 < pos2) {
                std::swap(a[pos1], a[pos2]);
            }
            increment(0);
        }
        return a;
    }

    BiFMIndexKStep(
        std::span<uint8_t const> _bwt,
        std::span<uint8_t const> _bwt_kstep,
        std::span<uint8_t const> _bwtRev_kstep,
        SparseArray _annotatedArray,
        VectorBool _annotatedArrayIsKStep)

        requires(!TReuseRev)
        : bwt_kstep{_bwt_kstep}
        , bwtRev_kstep{_bwtRev_kstep}
        , C{computeCSpan<Sigma>(_bwt)} //!TODO this should be possible with some smart call to bwt_kstep
        , C_kstep{computeC(bwtRev_kstep)}
        , CRev_kstep{computeC(bwt_kstep)}
        , annotatedArray{std::move(_annotatedArray)}
        , annotatedArrayIsKStep{std::move(_annotatedArrayIsKStep)}
    {
        if (
            _bwt.size() != bwtRev_kstep.size()
            || _bwt.size() != bwt_kstep.size()
        ) {
            throw std::runtime_error{"bwt don't have the same size: " + std::to_string(_bwt.size()) + " " + std::to_string(bwtRev_kstep.size()) + " " + std::to_string(bwt_kstep.size())};
        }
    }

    BiFMIndexKStep(std::span<uint8_t const> _bwt, std::span<uint8_t const> _bwt_kstep, SparseArray _annotatedArray, VectorBool _annotatedArrayIsKStep)
        requires(TReuseRev)
        : bwt_kstep{_bwt_kstep}
        , C{computeCSpan<Sigma>(_bwt)}
        , C_kstep{computeC(bwt_kstep)}
        , annotatedArray{std::move(_annotatedArray)}
        , annotatedArrayIsKStep{std::move(_annotatedArrayIsKStep)}
    {
    }


    /*
     * \param includeReversedInput assumes that the input data also has the reversed input data
     */
    BiFMIndexKStep(Sequence auto const& _sequence, SparseArray const& _annotatedSequence, auto _annotatedSequenceIsKStep, size_t _threadNbr, bool includeReversedInput = false) {
        if (_sequence.size() >= std::numeric_limits<size_t>::max()/2) {
            throw std::runtime_error{"sequence is longer than what this system is capable of handling"};
        }

        bool omegaSorting = !Delim_v; // Use omega sorting if no delimiter is being used

        // copy text into custom buffer
        auto inputText = createInputText(_sequence, omegaSorting, includeReversedInput);

        // create bwt, bwtRev and annotatedArray
//        auto [_bwt, _annotatedArray] = createBWTAndAnnotatedArray(inputText, _annotatedSequence, _threadNbr, omegaSorting);
        auto [_bwt, _bwt_kstep, _annotatedArray, _annotatedArrayIsKStep] = createBWTKStepAndAnnotatedArray<KStep, Sigma>(inputText, _annotatedSequence, _annotatedSequenceIsKStep, _threadNbr, omegaSorting);


        if constexpr (!TReuseRev) {
            #if defined(__GNUC__) && !defined(__clang__)
            #pragma GCC diagnostic push
            #pragma GCC diagnostic ignored "-Wstringop-overflow"
            std::ranges::reverse(inputText);
            #pragma GCC diagnostic pop

            #else
            std::ranges::reverse(inputText);

            #endif
            auto [_bwtRev, _bwtRev_kstep] = createBWTKStep<KStep, Sigma>(inputText, _threadNbr, omegaSorting);
            decltype(inputText){}.swap(inputText); // inputText memory can be deleted
            bwtRev_kstep = {_bwtRev_kstep};
        }

        // initialize this BiFMIndex properly
        bwt_kstep = {_bwt_kstep};
        C = computeCSpan<Sigma>(_bwt);//!TODO

        if constexpr (!TReuseRev) {
            C_kstep = computeC(bwtRev_kstep);
            CRev_kstep = computeC(bwt_kstep);

            // reorder CRev entries
            C_kstep = helperSwapC(C_kstep);
            CRev_kstep = helperSwapC(CRev_kstep);
        } else {
            C_kstep = computeC(bwt_kstep);
            C_kstep = helperSwapC(C_kstep);
        }

        annotatedArray        = std::move(_annotatedArray);
        annotatedArrayIsKStep = std::move(_annotatedArrayIsKStep);
    }


    /**!\brief Creates a BiFMIndex with a specified sampling rate
     *
     * \param _input a list of sequences
     * \param samplingRate rate of the sampling
     * \param includeReversedInput also adds all input and their reversed text
     */
    BiFMIndexKStep(Sequences auto const& _input, size_t samplingRate, size_t threadNbr, size_t seqOffset = 0, bool includeReversedInput = false) {
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

                    auto ret = std::optional<ADEntry>{};

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

                    auto ret = std::optional<ADEntry>{};

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
        auto annotatedSequenceIsKStep = VectorBool{};
        annotatedSequenceIsKStep.reserve(totalSize);
        for (size_t phase{0}; phase < totalSize; ++phase) {
            if (phase == 0) { // restarting
                refId = startRefId;
                pos = 0;
            }
            if (!includeReversedInput || phase*2 < totalSize) { // going forward
                assert(refId < inputSizes.size());
                assert(pos < inputSizes[refId]);

                auto lastSampledPosition = (pos / samplingRate) * samplingRate;
                auto distanceToLastSampledPosition = pos - lastSampledPosition;
                bool ret = (distanceToLastSampledPosition % KStep == 0);

                ++pos;
                if (inputSizes[refId] == pos) {
                    refId += 1;
                    pos = 0;
                }
                annotatedSequenceIsKStep.push_back(ret);
            } else { // going backwards
                assert(refId < inputSizes.size());
                assert(pos < inputSizes[refId]);

                auto lastSampledPosition = (pos / samplingRate) * samplingRate;
                auto distanceToLastSampledPosition = pos - lastSampledPosition;
                bool ret = (distanceToLastSampledPosition % KStep == 0);

                ++pos;
                if (inputSizes[refId] == pos) {
                    refId += 1;
                    pos = 0;
                }
                annotatedSequenceIsKStep.push_back(ret);
            }
        }
        *this = BiFMIndexKStep{inputText, annotatedSequence, annotatedSequenceIsKStep, threadNbr, /*includeReversedInput=*/false};
    }

    auto operator=(BiFMIndexKStep const&) -> BiFMIndexKStep& = delete;
    auto operator=(BiFMIndexKStep&& _other) noexcept -> BiFMIndexKStep& = default;

    size_t size() const {
        return bwt_kstep.size();
    }

    using LEntry = decltype(std::tuple_cat(DocumentEntries{}, std::tuple<size_t>{}));
    auto locate(size_t idx) const -> LEntry {
        uint64_t steps{};
        while (!annotatedArrayIsKStep[idx]) {
            auto symb = bwt_kstep.template symbol_limit<bitct>(idx);
            idx = bwt_kstep.template rank_limit<bitct>(idx, symb) + C[symb];
            steps += 1;
        }

        auto opt = annotatedArray.value(idx);
        while (!opt) {
            auto symb = bwt_kstep.symbol(idx);
            idx = bwt_kstep.rank(idx, symb) + C_kstep[symb];
            steps += KStep;
            opt = annotatedArray.value(idx);
        }
        return std::tuple_cat(*opt, std::tuple<size_t>{steps});
    }

    template <typename Archive>
    void serialize(this auto&& self, Archive& ar) {
        ar(/*self.bwt, */self.bwt_kstep, self.C, self.C_kstep, self.annotatedArray, self.annotatedArrayIsKStep);

        if constexpr (!std::same_as<RevBwtKStepType, std::nullptr_t>) {
            /*ar(self.bwtRev);*/
            ar(self.bwtRev_kstep);
            ar(self.CRev_kstep);
        }
    }
};

}
