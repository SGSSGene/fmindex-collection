#pragma once

#include "BitStack.h"
#include "Bitvector.h"
#include "cereal_tag.h"

#include <optional>

namespace fmindex_collection {

struct CSA {
    std::vector<uint64_t> ssa;
    Bitvector bv;
    size_t samplingRate; // distance between two samples (inside one sequence)

    CSA(std::vector<uint64_t> _ssa, BitStack const& bitstack, size_t _samplingRate)
        : ssa{std::move(_ssa)}
        , bv{bitstack.size, [&](size_t idx) {
            return bitstack.value(idx);
        }}
        , samplingRate{_samplingRate}
        {}
    CSA(CSA const&) = delete;
    CSA(CSA&& _other) noexcept = default;


    CSA(cereal_tag)
        : bv {cereal_tag{}}
    {}

    auto operator=(CSA const&) -> CSA& = delete;
    auto operator=(CSA&& _other) noexcept -> CSA& = default;

    size_t memoryUsage() const {
        return sizeof(ssa) + ssa.size() * sizeof(ssa.back())
            + bv.memoryUsage();
    }

    auto value(size_t idx) const -> std::optional<uint64_t> {
        if (!bv.value(idx)) {
            return std::nullopt;
        }
        return {ssa[bv.rank(idx)]};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(ssa, bv, samplingRate);
    }
};

}
#include "CSA32.h"
