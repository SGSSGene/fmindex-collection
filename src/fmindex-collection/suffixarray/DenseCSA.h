// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../BitStack.h"
#include "../bitvector/Bitvector.h"
#include "../bitvector/CompactBitvector.h"
#include "../DenseVector.h"
#include "concepts.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <tuple>

namespace fmindex_collection {

struct DenseCSA {
    DenseVector ssaPos;             // sampled suffix array - position inside the sequence
    DenseVector ssaSeq;             // sampled suffix array  - sequence id
    bitvector::CompactBitvector bv; // indicates if and which entry the sa → ssa
    size_t seqCount;

    static auto createJoinedCSA(DenseCSA const& lhs, DenseCSA const& rhs) -> DenseCSA {
        (void)lhs;
        (void)rhs;
        auto csa = DenseCSA{};
        csa.ssaPos = DenseVector(std::max(lhs.ssaPos.bits, rhs.ssaPos.bits));
        csa.ssaSeq = DenseVector(std::ceil(std::log2(lhs.seqCount + rhs.seqCount)));
        csa.seqCount = lhs.seqCount + rhs.seqCount;
        return csa;
    }


    DenseCSA() = default;
    DenseCSA(DenseCSA const&) = delete;
    DenseCSA(DenseCSA&&) noexcept = default;

    template <std::ranges::range Range>
        requires requires(Range r) {
            {*(r.begin())} -> std::same_as<std::optional<std::tuple<size_t, size_t>>>;
        }
    DenseCSA(Range _ssa, size_t sequencesCount, size_t longestSequence)
        : ssaPos(std::ceil(std::log2(longestSequence)))
        , ssaSeq(std::ceil(std::log2(sequencesCount)))
        , seqCount{sequencesCount}
    {
        for (auto o : _ssa) {
            bv.push_back(o.has_value());
            if (o) {
                auto [seqNr, pos] = *o;
                ssaPos.push_back(pos);
                ssaSeq.push_back(seqNr);
            }
        }
    }

    DenseCSA(DenseVector _ssaPos, DenseVector _ssaSeq, size_t _seqCount, BitStack const& bitstack)
        : ssaPos{std::move(_ssaPos)}
        , ssaSeq{std::move(_ssaSeq)}
        , bv{bitstack.size, [&](size_t idx) {
            return bitstack.value(idx);
        }}
        , seqCount{_seqCount} {
        assert(ssaPos.size() == ssaSeq.size());
    }

    DenseCSA(DenseVector _ssaPos, BitStack const& bitstack)
        : ssaPos{std::move(_ssaPos)}
        , ssaSeq(1)
        , bv{bitstack.size, [&](size_t idx) {
            return bitstack.value(idx);
        }}
        , seqCount{1} {
        for(size_t i{0}; i < ssaPos.size(); ++i) {
            ssaSeq.push_back(0);
        }
        assert(ssaPos.size() == ssaSeq.size());

    }

    DenseCSA(std::span<uint64_t const> sa, size_t samplingRate, std::span<size_t const> _inputSizes, bool reverse=false) {
        size_t bitsForSeqId = std::max(size_t{1}, size_t(std::ceil(std::log2(_inputSizes.size()))));
        assert(bitsForSeqId < 64);

        size_t largestText{};
        // Generate accumulated input
        auto accInputSizes = std::vector<uint64_t>{};
        accInputSizes.reserve(_inputSizes.size()+1);
        accInputSizes.emplace_back(0);
        for (auto len : _inputSizes) {
            accInputSizes.emplace_back(accInputSizes.back() + len);
            largestText = std::max(largestText, len);
        }

        // Construct sampled suffix array
        size_t bitsForPos   = std::max(size_t{1}, size_t(std::ceil(std::log2(largestText))));
        seqCount = _inputSizes.size();
        ssaPos = DenseVector(bitsForPos);
        ssaSeq = DenseVector(bitsForSeqId);
        ssaPos.reserve(sa.size() / samplingRate);
        ssaSeq.reserve(sa.size() / samplingRate);
        for (size_t i{0}; i < sa.size(); ++i) {
            bool sample = (sa[i] % samplingRate) == 0;
            if (sample) {
                // find subject id
                auto iter = std::upper_bound(accInputSizes.begin(), accInputSizes.end(), sa[i]);
                size_t subjId = std::distance(accInputSizes.begin(), iter) - 1;
                auto subjPos = sa[i] - accInputSizes[subjId];
                if (reverse) {
                    auto len = _inputSizes[subjId];
                    if (subjPos < len-1) {
                        subjPos = len - subjPos - 1;
                    } else {
                        subjPos = len;
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
    auto operator=(DenseCSA&&) noexcept -> DenseCSA& = default;

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

    void push_back(std::optional<std::tuple<size_t, size_t>> value) {
        bv.push_back(value.has_value());
        if (value) {
            auto [seqNr, pos] = *value;
            ssaPos.push_back(pos);
            ssaSeq.push_back(seqNr);
        }
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(ssaPos, ssaSeq, bv);
    }
};
static_assert(SuffixArray_c<DenseCSA>);


}
