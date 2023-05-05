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

    DenseCSA(std::span<int64_t const> sa, size_t _samplingRate, std::span<size_t const> _inputSizes, bool reverse=false)
        : ssaPos{cereal_tag{}}
        , ssaSeq{cereal_tag{}}
        , bv {cereal_tag{}}
        , samplingRate{_samplingRate}
    {
        size_t bitsForSeqId = std::max(size_t{1}, size_t(std::ceil(std::log2(_inputSizes.size()))));
        assert(bitsForSeqId < 64);
        size_t largestText  = *std::max_element(begin(_inputSizes), end(_inputSizes));
        size_t bitsForPos   = std::max(size_t{1}, size_t(std::ceil(std::log2(largestText))));

        // Generate accumulated input
        auto accInputSizes = std::vector<uint64_t>{};
        accInputSizes.reserve(_inputSizes.size()+1);
        accInputSizes.emplace_back(0);
        for (size_t i{0}; i < _inputSizes.size(); ++i) {
            accInputSizes.emplace_back(accInputSizes.back() + _inputSizes[i]);
        }

        // Annotate text with labels, naming the correct sequence id
        auto labels = std::vector<uint64_t>{};
        labels.reserve(sa.size() / samplingRate);

        for (size_t i{0}, subjId{0}; i < sa.size(); i += samplingRate) {
            while (i >= accInputSizes[subjId]) {
                subjId += 1;
            }
            labels.emplace_back(subjId-1);
        }

        // Construct sampled suffix array
        auto bitStack = fmindex_collection::BitStack{};
        ssaPos = DenseVector(bitsForPos);
        ssaSeq = DenseVector(bitsForSeqId);
        for (size_t i{0}; i < sa.size(); ++i) {
            bool sample = (sa[i] % samplingRate) == 0;
            bitStack.push(sample);
            if (sample) {
                auto subjId  = labels[sa[i] / samplingRate];
                auto subjPos = sa[i] - accInputSizes[subjId];
                if (reverse) {
                    subjPos = _inputSizes[subjId] - subjPos - 1;
                }
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
