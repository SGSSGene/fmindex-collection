// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "BitStack.h"
#include "Bitvector.h"
#include "cereal_tag.h"

#include <optional>

namespace fmindex_collection {

struct CSA32 {
    std::vector<uint32_t> ssa;
    Bitvector bv;

    CSA32(std::vector<uint64_t> const& _ssa, BitStack const& bitstack)
        : bv{bitstack.size, [&](size_t idx) {
            return bitstack.value(idx);
        }} {
            ssa.resize(_ssa.size());
            for (size_t i{0}; i < _ssa.size(); ++i) {
                ssa[i] = _ssa[i];
            }
        }
    CSA32(CSA32 const&) = delete;
    CSA32(CSA32&& _other) noexcept
        : ssa {std::move(_other.ssa)}
        , bv  {std::move(_other.bv)}
    {}


    CSA32(cereal_tag)
    {}

    auto operator=(CSA32 const&) -> CSA32& = delete;
    auto operator=(CSA32&& _other) noexcept -> CSA32& {
        ssa = std::move(_other.ssa);
        bv  = std::move(_other.bv);
        return *this;
    }


    size_t memoryUsage() const {
        return sizeof(ssa) + ssa.size() * sizeof(ssa.back())
            + bv.memoryUsage();
    }

    auto value(size_t idx) const -> std::optional<uint64_t> {
        if (!bv.symbol(idx)) {
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
