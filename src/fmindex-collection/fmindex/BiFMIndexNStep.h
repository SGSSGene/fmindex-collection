// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
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


template <size_t TSigma, template <size_t> typename String = string::FlattenedBitvectors_512_64k, SparseArray_c SparseArray = suffixarray::SparseArray<std::tuple<uint32_t, uint32_t>>, bool TDelim=true, bool TReuseRev=false, size_t TNStep=2>
    requires String_c<String<TSigma>>
struct BiFMIndexNStep {
    using ADEntry = SparseArray::value_t;

    using NoDelim = BiFMIndexNStep<TSigma, String, SparseArray, false, TReuseRev, TNStep>;
    using ReuseRev = BiFMIndexNStep<TSigma, String, SparseArray, TDelim, true, TNStep>;

    template <size_t TNewNStep>
    using SetNStep = BiFMIndexNStep<TSigma, String, SparseArray, TDelim, TReuseRev, TNewNStep>;

    static size_t constexpr Sigma     = TSigma;
    static size_t constexpr FirstSymb = TDelim?1:0;
    static bool constexpr Delim_v     = TDelim;
    static bool constexpr ReuseRev_v  = TReuseRev;


    static size_t constexpr NStep = TNStep;
    static size_t constexpr SigmaNStep = my_pow_base<Sigma>(NStep);

    // Set RevBwtType to std::nullptr_t to indicate that it should not be used
    using RevBwtType = std::conditional_t<TReuseRev, std::nullptr_t, String<Sigma>>;
    using RevBwtNStepType = std::conditional_t<TReuseRev, std::nullptr_t, String<SigmaNStep>>;
    static_assert(
        std::same_as<RevBwtType, std::nullptr_t> == std::same_as<RevBwtNStepType, std::nullptr_t>
        , "either RevBWTType and RevBwtNStepType have nullptr_t or none"
    );


    String<Sigma> bwt;
    RevBwtType bwtRev;
    String<SigmaNStep> bwt_nstep;
    RevBwtNStepType bwtRev_nstep;

    std::array<size_t, Sigma+1> C{};
    std::array<size_t, SigmaNStep+1> C_nstep{};
    std::array<size_t, SigmaNStep+1> CRev_nstep{};
    SparseArray annotatedArray;

    BiFMIndexNStep() = default;
    BiFMIndexNStep(BiFMIndexNStep&&) noexcept = default;


    static auto helperSwapC(std::array<size_t, SigmaNStep+1> a) {
        auto values = std::array<size_t, NStep>{};
        auto increment = std::function<void(size_t pos)>{};

        increment = [&](size_t pos) {
            if (pos == NStep) return;
            values[pos] += 1;
            if (values[pos] == Sigma) {
                values[pos] = 0;
                increment(pos+1);
            }
        };
        auto asValueFwd = [](std::array<size_t, NStep> const& v) {
            size_t value{};
            for (size_t i{0}; i < NStep; ++i) {
                value = value * Sigma + v[i];
            }
            return value;
        };
        auto asValueBwd = [](std::array<size_t, NStep> const& v) {
            size_t value{};
            for (size_t i{0}; i < NStep; ++i) {
                value = value * Sigma + v[NStep -1 -i];
            }
            return value;
        };

        for (size_t i{0}; i < std::pow(Sigma, NStep); ++i) {
            auto pos1 = asValueFwd(values);
            auto pos2 = asValueBwd(values);
            if (pos1 < pos2) {
                std::swap(a[pos1], a[pos2]);
            }
            increment(0);
        }
        return a;
    }

    BiFMIndexNStep(
        std::span<uint8_t const> _bwt,
        std::span<uint8_t const> _bwtRev,
        std::span<uint8_t const> _bwt_nstep,
        std::span<uint8_t const> _bwtRev_nstep,
        SparseArray _annotatedArray)

        requires(!TReuseRev)
        : bwt{_bwt}
        , bwtRev{_bwtRev}
        , bwt_nstep{_bwt_nstep}
        , bwtRev_nstep{_bwtRev_nstep}
        , C{computeC(bwt)}
        , C_nstep{computeC(bwtRev_nstep)}
        , CRev_nstep{computeC(bwt_nstep)}
        , annotatedArray{std::move(_annotatedArray)}
    {
        assert(bwt.size() == bwtRev.size());
        if (
            bwt.size() != bwtRev.size()
            || bwt.size() != bwt_nstep.size()
            || bwt.size() != bwtRev_nstep.size()
        ) {
            throw std::runtime_error("bwt don't have the same size: " + std::to_string(bwt.size()) + " " + std::to_string(bwtRev.size()));
        }
    }

    BiFMIndexNStep(std::span<uint8_t const> _bwt, std::span<uint8_t const> _bwt_nstep, SparseArray _annotatedArray)
        requires(TReuseRev)
        : bwt{_bwt}
        , bwt_nstep{_bwt_nstep}
        , C{computeC(bwt)}
        , C_nstep{computeC(bwt_nstep)}
        , annotatedArray{std::move(_annotatedArray)}
    {}


    /*
     * \param includeReversedInput assumes that the input data also has the reversed input data
     */
    BiFMIndexNStep(Sequence auto const& _sequence, SparseArray const& _annotatedSequence, size_t _threadNbr, bool includeReversedInput = false) {
        if (_sequence.size() >= std::numeric_limits<size_t>::max()/2) {
            throw std::runtime_error{"sequence is longer than what this system is capable of handling"};
        }

        bool omegaSorting = !Delim_v; // Use omega sorting if no delimiter is being used

        // copy text into custom buffer
        auto inputText = createInputText(_sequence, omegaSorting, includeReversedInput);

        // create bwt, bwtRev and annotatedArray
//        auto [_bwt, _annotatedArray] = createBWTAndAnnotatedArray(inputText, _annotatedSequence, _threadNbr, omegaSorting);
        auto [_bwt, _bwt_nstep, _annotatedArray] = createBWTNStepAndAnnotatedArray<NStep, Sigma>(inputText, _annotatedSequence, _threadNbr, omegaSorting);


        if constexpr (!TReuseRev) {
            #if defined(__GNUC__) && !defined(__clang__)
            #pragma GCC diagnostic push
            #pragma GCC diagnostic ignored "-Wstringop-overflow"
            std::ranges::reverse(inputText);
            #pragma GCC diagnostic pop

            #else
            std::ranges::reverse(inputText);

            #endif
            auto [_bwtRev, _bwtRev_nstep] = createBWTNStep<NStep, Sigma>(inputText, _threadNbr, omegaSorting);
            decltype(inputText){}.swap(inputText); // inputText memory can be deleted
            bwtRev       = {_bwtRev};
            bwtRev_nstep = {_bwtRev_nstep};
        }

        // initialize this BiFMIndex properly
        bwt = {_bwt};
        bwt_nstep = {_bwt_nstep};
        C = computeC(bwt);

        if constexpr (!TReuseRev) {
            C_nstep = computeC(bwtRev_nstep);
            CRev_nstep = computeC(bwt_nstep);

            // reorder CRev entries
            C_nstep = helperSwapC(C_nstep);
            CRev_nstep = helperSwapC(CRev_nstep);
        } else {
            C_nstep = computeC(bwt_nstep);
            C_nstep = helperSwapC(C_nstep);
        }

        annotatedArray = std::move(_annotatedArray);
    }


    /**!\brief Creates a BiFMIndex with a specified sampling rate
     *
     * \param _input a list of sequences
     * \param samplingRate rate of the sampling
     * \param includeReversedInput also adds all input and their reversed text
     */
    BiFMIndexNStep(Sequences auto const& _input, size_t samplingRate, size_t threadNbr, size_t seqOffset = 0, bool includeReversedInput = false) {
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

        *this = BiFMIndexNStep{inputText, annotatedSequence, threadNbr, /*includeReversedInput=*/false};
    }

    auto operator=(BiFMIndexNStep const&) -> BiFMIndexNStep& = delete;
    auto operator=(BiFMIndexNStep&& _other) noexcept -> BiFMIndexNStep& = default;

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
        ar(self.bwt, self.bwt_nstep, self.C, self.C_nstep, self.annotatedArray);

        if constexpr (!std::same_as<RevBwtType, std::nullptr_t>) {
            ar(self.bwtRev);
            ar(self.bwtRev_nstep);
            ar(self.CRev_nstep);
        }
    }
};

}
