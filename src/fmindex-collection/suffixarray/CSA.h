// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once


#include "../bitvector/L0L1_NBitvector.h"
#include "concepts.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <optional>
#include <tuple>


namespace fmindex_collection {

namespace suffixarray {
template <typename T, typename SAEntry, Bitvector_c Bitvector>
auto createSampling(std::vector<SAEntry> const& sa,
                    Bitvector const& textAnnotationValid,
                    std::vector<T> const& textAnnotation) -> std::pair<bitvector::L0L1_512_64kBitvector, std::vector<T>> {
    auto bv  = bitvector::L0L1_512_64kBitvector{};
    auto ssa = std::vector<T>{};
    for (size_t i{0}; i < sa.size(); ++i) {
        auto textPos = sa[i];
        auto valid = textAnnotationValid.symbol(textPos);
        bv.push_back(valid);
        if (valid) {
            auto rank = textAnnotationValid.rank(textPos);
            ssa.push_back(textAnnotation[rank]);
        }
    }
    return {bv, ssa};
}
}

struct CSA {
    using Bitvector = bitvector::L0L1_512_64kBitvector;
    std::vector<uint64_t> ssa;
    Bitvector             bv;
    size_t                bitsForPosition;   // bits reserved for position
    uint64_t              bitPositionMask;   // Bit mask, to extract the position from ssa
    size_t                seqCount;          // Number of sequences


    static auto createJoinedCSA(CSA const& lhs, CSA const& rhs) -> CSA {
        auto csa = CSA{};
        size_t bitsForPosition = std::max(lhs.bitsForPosition, rhs.bitsForPosition);
        auto seqBits = size_t(std::ceil(std::log2(lhs.seqCount + rhs.seqCount)));
        if (seqBits + bitsForPosition > 64) {
            throw std::runtime_error{"Can't merge indices, requires more than 64bit"};
        }

        csa.bitsForPosition = bitsForPosition;
        csa.bitPositionMask = (1ull<<bitsForPosition)-1;
        return csa;
    }


    CSA() = default;
    CSA(CSA const&) = delete;
    CSA(CSA&&) noexcept = default;

    template <std::ranges::range Range>
        requires requires(Range r) {
            {*(r.begin())} -> std::same_as<std::optional<std::tuple<size_t, size_t>>>;
        }
    CSA(Range _ssa, size_t sequencesCount, size_t longestSequence)
        : seqCount{sequencesCount}
    {
        bitsForPosition = size_t(std::ceil(std::log2(longestSequence)));
        size_t bitsForSeqId = longestSequence;
        if (bitsForPosition + bitsForSeqId > 64) {
            throw std::runtime_error{"requires more than 64bit to encode sequence length and number of sequence"};
        }
        bitPositionMask = (uint64_t{1}<<bitsForPosition)-1;

        for (auto o : _ssa) {
            push_back(o);
        }
    }

    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint8_t>
    CSA(std::vector<uint64_t> _ssa, range_t const& bitstack, size_t _bitsForPosition, size_t _seqCount)
        : ssa{std::move(_ssa)}
        , bv{bitstack}
        , bitsForPosition{_bitsForPosition}
        , bitPositionMask{(uint64_t{1}<<bitsForPosition)-1}
        , seqCount{_seqCount}
    {}

    template <typename T>
    CSA(std::vector<T> const& sa, size_t samplingRate, std::span<size_t const> _inputSizes, bool reverse=false, size_t seqOffset=0)
        : seqCount{_inputSizes.size()}
    {
        assert(samplingRate != 0);
        size_t longestSequence = std::accumulate(_inputSizes.begin(), _inputSizes.end(), size_t{}, [](size_t lhs, auto rhs) {
            return std::max(lhs, rhs);
        });
        bitsForPosition = size_t(std::ceil(std::log2(longestSequence)));
        size_t bitsForSeqId = std::max(size_t{1}, size_t(std::ceil(std::log2(_inputSizes.size() + seqOffset))));
        if (bitsForPosition + bitsForSeqId > 64) {
            throw std::runtime_error{"requires more than 64bit to encode sequence length and number of sequence"};
        }
        bitPositionMask = (uint64_t{1}<<bitsForPosition)-1;


        // create a sampling in text space
        auto textAnnotationValid = Bitvector{};
        auto textAnnotation = std::vector<uint64_t>{};
        for (size_t refId{0}; refId < _inputSizes.size(); ++refId) {
            for (size_t posId{0}; posId < _inputSizes[refId]; ++posId) {
                auto sample = (posId % samplingRate) == 0;
                textAnnotationValid.push_back(sample);
                if (!sample) continue;
                if (!reverse) {
                    textAnnotation.emplace_back(((refId+seqOffset) << bitsForPosition) | posId);
                } else {
                    auto len = _inputSizes[refId];
                    textAnnotation.emplace_back(((refId+seqOffset) << bitsForPosition) | (len - posId - 1));
                }
            }
        }
        // converts text space into sampled suffix array space
        std::tie(bv, ssa) = suffixarray::createSampling(sa, textAnnotationValid, textAnnotation);
    }

    auto operator=(CSA const&) -> CSA& = delete;
    auto operator=(CSA&&) noexcept -> CSA& = default;

    size_t memoryUsage() const {
        return sizeof(ssa) + ssa.size() * sizeof(ssa.back());
    }

    auto value(size_t idx) const -> std::optional<std::tuple<uint64_t, uint64_t>> {
        if (!bv.symbol(idx)) {
            return std::nullopt;
        }
        auto v = ssa[bv.rank(idx)];
        auto chr = v >> bitsForPosition;
        auto pos = v & bitPositionMask;

        return std::make_tuple(chr, pos);
    }

    void push_back(std::optional<std::tuple<size_t, size_t>> value) {
        bv.push_back(value.has_value());
        if (value) {
            auto [seqNr, pos] = *value;
            ssa.push_back((seqNr << bitsForPosition) + pos);
        }
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(ssa, bv, bitsForPosition, bitPositionMask, seqCount);
    }
};
static_assert(SuffixArray_c<CSA>);

}
