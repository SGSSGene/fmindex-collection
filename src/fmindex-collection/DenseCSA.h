#pragma once

#include "BitStack.h"
#include "Bitvector.h"
#include "BitvectorCompact.h"
#include "DenseVector.h"
#include "cereal_tag.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <tuple>

namespace fmindex_collection {

struct DenseCSA {
    DenseVector ssaPos;     // sampled suffix array - position inside the sequence
    DenseVector ssaSeq;     // sampled suffix array  - sequence id
    BitvectorCompact bv;    // indicates if and which entry the sa → ssa
    size_t samplingRate;    // distance between two samples (inside one sequence)

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
        , bv {cereal_tag{}}
    {}

    DenseCSA(std::vector<int64_t> const& sa, size_t _samplingRate, std::vector<size_t> const& _inputSizes, bool reverse=false)
        : ssaPos{cereal_tag{}}
        , ssaSeq{cereal_tag{}}
        , bv {cereal_tag{}}
        , samplingRate{_samplingRate}
    {
        size_t bitsForSeqId = std::max(1ul, size_t(std::ceil(std::log2(_inputSizes.size()))));
        assert(bitsForSeqId < 64);
        size_t largestText  = *std::max_element(begin(_inputSizes), end(_inputSizes));
        size_t bitsForPos   = std::max(1ul, size_t(std::ceil(std::log2(largestText))));


        auto bitStack = fmindex_collection::BitStack{};
        auto ssa      = std::vector<uint64_t>{};
        if (samplingRate > 0) {
            ssa.reserve(sa.size() / samplingRate + 1 + _inputSizes.size());
        }

        ssaPos = DenseVector(bitsForPos);
        ssaSeq = DenseVector(bitsForSeqId);

        auto accIter = _inputSizes.begin();
        size_t subjId{};
        size_t subjPos{};

        size_t lastSamplingPos{};

        auto newLabels = std::vector<std::tuple<uint64_t, uint64_t>>{};
        newLabels.resize(sa.size(), std::make_tuple(std::numeric_limits<uint64_t>::max(), 0ul));

        for (size_t i{0}; i < newLabels.size(); ++i, ++subjPos) {
            while (subjPos >= *accIter) {
                subjPos -= *accIter;
                ++subjId;
                ++accIter;
            }
            bool sample = (samplingRate == 0) || ((i-lastSamplingPos) % samplingRate == 0) || (subjPos == 0);
            if (sample) {
                lastSamplingPos = i;
                auto pos = subjPos;
                if (reverse) {
                    pos = _inputSizes[subjId] - pos -1;
                }
                newLabels[i] = {subjId, pos};
            }
        }

        for (size_t i{0}; i < sa.size(); ++i) {
            auto [subjId, subjPos] = newLabels[sa[i]];
            bool sample = subjId != std::numeric_limits<uint64_t>::max();
            bitStack.push(sample);
            if (sample) {
                ssaSeq.push_back(subjId);
                ssaPos.push_back(subjPos);
            }
        }


        this->bv  = BitvectorCompact{bitStack.size, [&](size_t idx) {
            return bitStack.value(idx);
        }};
    }


    auto operator=(DenseCSA const&) -> DenseCSA& = delete;
    auto operator=(DenseCSA&& _other) noexcept -> DenseCSA& = default;

    size_t memoryUsage() const {
        return (ssaPos.bit_size() + ssaSeq.bit_size()) / 8 + bv.memoryUsage();
    }

    auto value(size_t idx) const -> std::optional<std::tuple<uint64_t, uint64_t>> {
        if (!bv.value(idx)) {
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

}
