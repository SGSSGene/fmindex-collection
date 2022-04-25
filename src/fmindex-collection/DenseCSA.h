#pragma once

#include "BitStack.h"
#include "Bitvector.h"
#include "BitvectorCompact.h"
#include "DenseVector.h"
#include "cereal_tag.h"

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
