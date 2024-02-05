// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "BitStack.h"
#include "DenseVector.h"
#include "bitvector/Bitvector.h"
#include "bitvector/CompactBitvector.h"
#include "cereal_tag.h"
#include "suffixarray/concepts.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <tuple>

namespace fmindex_collection {

struct DenseCSA {
    DenseVector ssaPos;             // sampled suffix array - position inside the sequence
    DenseVector ssaSeq;             // sampled suffix array  - sequence id
    bitvector::CompactBitvector bv; // indicates if and which entry the sa → ssa
    size_t samplingRate;    // distance between two samples (inside one sequence)

    DenseCSA() = default;
    DenseCSA(DenseVector _ssaPos, DenseVector _ssaSeq, BitStack const& bitstack, size_t _samplingRate)
        : ssaPos{std::move(_ssaPos)}
        , ssaSeq{std::move(_ssaSeq)}
        , bv{bitstack.size, [&](size_t idx) {
            return bitstack.value(idx);
        }}
        , samplingRate{_samplingRate}
        {
        assert(ssaPos.size() == ssaSeq.size());
    }

    DenseCSA(DenseVector _ssaPos, BitStack const& bitstack, size_t _samplingRate)
        : ssaPos{std::move(_ssaPos)}
        , ssaSeq(1)
        , bv{bitstack.size, [&](size_t idx) {
            return bitstack.value(idx);
        }}
        , samplingRate{_samplingRate}
        {
        for(size_t i{0}; i < ssaPos.size(); ++i) {
            ssaSeq.push_back(0);
        }
        assert(ssaPos.size() == ssaSeq.size());

    }

    DenseCSA(DenseCSA const&) = delete;
    DenseCSA(DenseCSA&& _other) noexcept = default;


    DenseCSA(cereal_tag)
        : ssaPos{cereal_tag{}}
        , ssaSeq{cereal_tag{}}
    {}

    DenseCSA(std::span<uint64_t const> sa, size_t _samplingRate, std::span<std::tuple<size_t, size_t> const> _inputSizes, bool reverse=false)
        : ssaPos{cereal_tag{}}
        , ssaSeq{cereal_tag{}}
        , samplingRate{_samplingRate}
    {
        size_t bitsForSeqId = std::max(size_t{1}, size_t(std::ceil(std::log2(_inputSizes.size()))));
        assert(bitsForSeqId < 64);

        size_t largestText{};
        // Generate accumulated input
        auto accInputSizes = std::vector<uint64_t>{};
        accInputSizes.reserve(_inputSizes.size()+1);
        accInputSizes.emplace_back(0);
        for (size_t i{0}; i < _inputSizes.size(); ++i) {
            auto [len, delCt] = _inputSizes[i];
            accInputSizes.emplace_back(accInputSizes.back() + len + delCt);
            largestText = std::max(largestText, len);
        }

        // Construct sampled suffix array
        size_t bitsForPos   = std::max(size_t{1}, size_t(std::ceil(std::log2(largestText))));

        ssaPos = DenseVector(bitsForPos);
        ssaSeq = DenseVector(bitsForSeqId);
        ssaPos.reserve(sa.size() / _samplingRate);
        ssaSeq.reserve(sa.size() / _samplingRate);
        for (size_t i{0}; i < sa.size(); ++i) {
            bool sample = (sa[i] % samplingRate) == 0;
            if (sample) {
                // find subject id
                auto iter = std::upper_bound(accInputSizes.begin(), accInputSizes.end(), sa[i]);
                size_t subjId = std::distance(accInputSizes.begin(), iter) - 1;
                auto subjPos = sa[i] - accInputSizes[subjId];
                if (reverse) {
                    auto [len, delCt] = _inputSizes[subjId];
                    if (subjPos < len) {
                        subjPos = len - subjPos;
                    } else {
                        subjPos = len+1;
                    }
                }
                ssaSeq.push_back(subjId);
                ssaPos.push_back(subjPos);
            }
        }

        this->bv  = bitvector::CompactBitvector{sa.size(), [&](size_t idx) {
            return (sa[idx] % samplingRate) == 0;
        }};
    }


    auto operator=(DenseCSA const&) -> DenseCSA& = delete;
    auto operator=(DenseCSA&& _other) noexcept -> DenseCSA& = default;

    size_t memoryUsage() const {
        return (ssaPos.bit_size() + ssaSeq.bit_size()) / 8;
    }

    auto value(size_t idx) const -> std::optional<std::tuple<uint64_t, uint64_t>> {
        if (!bv.symbol(idx)) {
            return std::nullopt;
        }
        auto rank = bv.rank(idx);
        return std::make_tuple(ssaSeq[rank], ssaPos[rank]);
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(ssaPos, ssaSeq, bv, samplingRate);
    }
};
static_assert(SuffixArray_c<DenseCSA>);


}
