#pragma once

#include "BitStack.h"
#include "Bitvector.h"

#include <optional>

struct CSA {
    std::vector<uint64_t> ssa;
    Bitvector bv;

    CSA(std::vector<uint64_t> _ssa, BitStack const& bitstack)
        : ssa{std::move(_ssa)}
        , bv{bitstack.size, [&](size_t idx) {
            return bitstack.value(idx);
        }}
    {}

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

};
