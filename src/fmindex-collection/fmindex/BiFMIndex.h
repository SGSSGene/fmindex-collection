// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../string/concepts.h"
#include "../suffixarray/CSA.h"
#include "../suffixarray/SparseArray.h"
#include "../utils.h"

#include <algorithm>
#include <cassert>

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

template <String_c String, SuffixArray_c TCSA = CSA, typename T = std::tuple<size_t, size_t>>
struct BiFMIndex {
    using ADEntry = T;

    static size_t constexpr Sigma = String::Sigma;

    String bwt;
    String bwtRev;
    std::array<size_t, Sigma+1> C{};
    suffixarray::SparseArray<T> annotatedArray;

    BiFMIndex() = default;
    BiFMIndex(BiFMIndex&&) noexcept = default;

    BiFMIndex(std::span<uint8_t const> _bwt, std::span<uint8_t const> _bwtRev, TCSA _csa)
        : bwt{_bwt}
        , bwtRev{_bwtRev}
        , annotatedArray{std::views::iota(size_t{0}, _csa.bv.size()) | std::views::transform([&](size_t i) {
            return _csa.value(i);
        })}
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

    BiFMIndex(std::span<uint8_t const> _bwt, std::span<uint8_t const> _bwtRev, suffixarray::SparseArray<T> _annotatedArray)
        : bwt{_bwt}
        , bwtRev{_bwtRev}
        , annotatedArray{std::move(_annotatedArray)}
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

    BiFMIndex(Sequence auto const& _sequence, suffixarray::SparseArray<T> const& _annotatedSequence, size_t _threadNbr, bool omegaSorting = false) {
        if (_sequence.size() >= std::numeric_limits<size_t>::max()/2) {
            throw std::runtime_error{"sequence is longer than what this system is capable of handling"};
        }

        // copy text into custom buffer
        auto inputText = std::vector<uint8_t>{};
        auto inputTextSpan = std::span<uint8_t const>{};
        if (omegaSorting) {
            inputText.resize(_sequence.size()*2);
            for (size_t i{0}; i < _sequence.size(); ++i) {
                inputText[i] = _sequence[i];
                inputText[i+_sequence.size()] = _sequence[i];
            }
            inputTextSpan = {inputText.begin(), inputText.begin() + inputText.size()/2};
        } else {
            inputText.resize(_sequence.size());
            for (size_t i{0}; i < _sequence.size(); ++i) {
                inputText[i] = _sequence[i];
            }
            inputTextSpan = {inputText.begin(), inputText.begin() + inputText.size()};
        }

        auto removeInvalidSA = [&](auto& sa) {
            if (omegaSorting) { // using omega sorting, remove half of the entries
                auto halfSize = inputText.size() / 2;
                auto [first, last] = std::ranges::remove_if(sa, [&](auto e) {
                    return e >= halfSize;
                });
                sa.erase(first, last);
            }
        };

        if (inputText.size() < std::numeric_limits<int32_t>::max()) { // only 32bit SA required
            // create BurrowsWheelerTransform and CompressedSuffixArray
            auto [bwt, annotatedArray] = [&]() {
                auto sa  = createSA32(inputText, _threadNbr);
                removeInvalidSA(sa);

                auto bwt = createBWT32(inputTextSpan, sa);
                auto annotatedArray = suffixarray::SparseArray<T>{
                    sa | std::views::transform([&](size_t i) -> std::optional<T> {
                        return _annotatedSequence.value(i);
                    })
                };
                return std::make_tuple(std::move(bwt), std::move(annotatedArray));
            }();

            // create BurrowsWheelerTransform on reversed text
            auto bwtRev = [&]() {
                //std::ranges::reverse(inputText); //!TODO produces a weird error
                for (size_t i{0}; i < inputText.size()/2; ++i) {
                    std::swap(inputText[i], inputText[inputText.size() - i - 1]);
                }
                auto saRev  = createSA32(inputText, _threadNbr);
                removeInvalidSA(saRev);

                auto bwtRev = createBWT32(inputTextSpan, saRev);
                return bwtRev;
            }();

            decltype(inputText){}.swap(inputText); // inputText memory can be deleted

            *this = BiFMIndex{bwt, bwtRev, std::move(annotatedArray)};
        } else { // required 64bit SA required
            // create BurrowsWheelerTransform and CompressedSuffixArray
            auto [bwt, annotatedArray] = [&]() {
                auto sa  = createSA64(inputText, _threadNbr);
                removeInvalidSA(sa);

                auto bwt = createBWT64(inputTextSpan, sa);
                auto annotatedArray = suffixarray::SparseArray<T>{
                    sa | std::views::transform([&](size_t i) -> std::optional<T> {
                        return _annotatedSequence.value(i);
                    })
                };
                return std::make_tuple(std::move(bwt), std::move(annotatedArray));
            }();

            // create BurrowsWheelerTransform on reversed text
            auto bwtRev = [&]() {
                //std::ranges::reverse(inputText); //!TODO produces a weird error
                for (size_t i{0}; i < inputText.size()/2; ++i) {
                    std::swap(inputText[i], inputText[inputText.size() - i - 1]);
                }
                auto saRev  = createSA64(inputText, _threadNbr);
                removeInvalidSA(saRev);
                auto bwtRev = createBWT64(inputTextSpan, saRev);
                return bwtRev;
            }();
            (void)bwt;
            (void)annotatedArray;
            (void)bwtRev;

            decltype(inputText){}.swap(inputText); // inputText memory can be deleted

            *this = BiFMIndex{bwt, bwtRev, std::move(annotatedArray)};
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

        size_t refId{0};
        size_t pos{0};

        //!TODO what about empty strings?
        while (inputSizes[refId] == pos) {
            refId += 1;
            assert(inputSizes.size() < refId);
        }

        auto annotatedSequence = suffixarray::SparseArray<T>{
            std::views::iota(size_t{0}, totalSize) | std::views::transform([&](size_t) -> std::optional<T> {
                assert(refId < inputSizes.size());
                assert(pos < inputSizes[refId]);

                auto ret = std::optional<T>{std::nullopt};

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

    auto locate(size_t idx) const -> std::tuple<T, size_t> {
        if constexpr (requires(String t) {{ t.hasValue(size_t{}) }; }) {
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
                if constexpr (requires(String t) { { t.rank_symbol(size_t{}) }; }) {
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

    auto single_locate_step(size_t idx) const -> std::optional<T> {
        return annotatedArray.value(idx);
    }


    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bwt, bwtRev, C, annotatedArray);
    }
};

}
#ifdef WORKAROUND_GCC_MACOS14
    #pragma GCC diagnostic pop
#endif
