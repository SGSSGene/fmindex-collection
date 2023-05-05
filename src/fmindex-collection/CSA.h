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
        , bitPositionMask{(1ull<<bitsForPosition)-1}
    {}
    CSA(CSA const&) = delete;
    CSA(CSA&& _other) noexcept = default;

    CSA(cereal_tag)
        : bv {cereal_tag{}}
    {}

    CSA(std::span<int64_t const> sa, size_t _samplingRate, std::span<size_t const> _inputSizes, bool reverse=false)
        : samplingRate{_samplingRate}
    {
        size_t bitsForSeqId = std::max(size_t{1}, size_t(std::ceil(std::log2(_inputSizes.size()))));
        assert(bitsForSeqId < 64);

        bitsForPosition = 64 - bitsForSeqId;
        bitPositionMask = (1ull<<bitsForPosition)-1;

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
        auto ssa = std::vector<uint64_t>{};
        ssa.reserve(_inputSizes.size() / _samplingRate);
        for (size_t i{0}; i < sa.size(); ++i) {
            bool sample = (sa[i] % samplingRate) == 0;
            bitStack.push(sample);
            if (sample) {
                auto subjId  = labels[sa[i] / samplingRate];
                auto subjPos = sa[i] - accInputSizes[subjId];
                if (reverse) {
                    subjPos = _inputSizes[subjId] - subjPos - 1;
                }
                ssa.emplace_back(subjPos | (subjId << bitsForPosition));
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
