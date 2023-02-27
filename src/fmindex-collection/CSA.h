#pragma once

#include "BitStack.h"
#include "Bitvector.h"
#include "BitvectorCompact.h"
#include "cereal_tag.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <optional>
#include <tuple>


namespace fmindex_collection {

struct CSA {
    std::vector<uint64_t> ssa;
    BitvectorCompact bv;
    size_t samplingRate;    // distance between two samples (inside one sequence)
    size_t bitsForPosition; // bits reserved for position
    size_t bitPositionMask;


    CSA(std::vector<uint64_t> _ssa, BitStack const& bitstack, size_t _samplingRate, size_t _bitsForPosition)
        : ssa{std::move(_ssa)}
        , bv{bitstack.size, [&](size_t idx) {
            return bitstack.value(idx);
        }}
        , samplingRate{_samplingRate}
        , bitsForPosition{_bitsForPosition}
        , bitPositionMask{(1ul<<bitsForPosition)-1}
    {}
    CSA(CSA const&) = delete;
    CSA(CSA&& _other) noexcept = default;

    CSA(cereal_tag)
        : bv {cereal_tag{}}
    {}

    CSA(std::vector<int64_t> const& sa, size_t _samplingRate, std::vector<size_t> const& _inputSizes, bool reverse=false)
        : samplingRate{_samplingRate}
    {
        size_t bitsForSeqId = std::max(size_t{1}, size_t(std::ceil(std::log2(_inputSizes.size()))));
        assert(bitsForSeqId < 64);

        bitsForPosition = 64 - bitsForSeqId;
        bitPositionMask = (1ul<<bitsForPosition)-1;


        auto bitStack = fmindex_collection::BitStack{};
        auto ssa      = std::vector<uint64_t>{};
        if (samplingRate > 0) {
            ssa.reserve(sa.size() / samplingRate + 1 + _inputSizes.size());
        }

        auto newLabels = std::vector<uint64_t>{};
        newLabels.resize(sa.size(), std::numeric_limits<uint64_t>::max());

        auto accIter = _inputSizes.begin();
        size_t subjId{};
        size_t subjPos{};

        size_t lastSamplingPos{};
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
                newLabels[i] = pos | (subjId << bitsForPosition);
            }
        }

        for (size_t i{0}; i < sa.size(); ++i) {
            bool sample = newLabels[sa[i]] != std::numeric_limits<uint64_t>::max();
            bitStack.push(sample);
            if (sample) {
                ssa.push_back(newLabels[sa[i]]);
            }
        }
        this->ssa = std::move(ssa);
        this->bv  = BitvectorCompact{bitStack.size, [&](size_t idx) {
            return bitStack.value(idx);
        }};
    }


    auto operator=(CSA const&) -> CSA& = delete;
    auto operator=(CSA&& _other) noexcept -> CSA& = default;

    size_t memoryUsage() const {
        return sizeof(ssa) + ssa.size() * sizeof(ssa.back())
            + bv.memoryUsage();
    }

    auto value(size_t idx) const -> std::optional<std::tuple<uint64_t, uint64_t>> {
        if (!bv.value(idx)) {
            return std::nullopt;
        }
        auto v = ssa[bv.rank(idx)];
        auto chr = v >> bitsForPosition;
        auto pos = v & bitPositionMask;

        return std::make_tuple(chr, pos);
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(ssa, bv, samplingRate, bitsForPosition, bitPositionMask);
    }
};

}
