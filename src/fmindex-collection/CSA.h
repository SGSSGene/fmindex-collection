#pragma once

#include "BitStack.h"
#include "Bitvector.h"
#include "cereal_tag.h"

#include <optional>

namespace fmindex_collection {

struct CSA {
    std::vector<uint64_t> ssa;
    Bitvector bv;

    CSA(std::vector<uint64_t> _ssa, BitStack const& bitstack)
        : ssa{std::move(_ssa)}
        , bv{bitstack.size, [&](size_t idx) {
            return bitstack.value(idx);
        }} {}
    CSA(CSA const&) = delete;
    CSA(CSA&& _other) noexcept
        : ssa {std::move(_other.ssa)}
        , bv  {std::move(_other.bv)}
    {}


    CSA(cereal_tag)
        : bv {cereal_tag{}}
    {}

    auto operator=(CSA const&) -> CSA& = delete;
    auto operator=(CSA&& _other) noexcept -> CSA& {
        ssa = std::move(_other.ssa);
        bv  = std::move(_other.bv);
        return *this;
    }


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
        ar(ssa, bv);
    }
};

}
