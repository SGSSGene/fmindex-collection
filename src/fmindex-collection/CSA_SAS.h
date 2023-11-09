// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "cereal_tag.h"

#include <cstdint>
#include <optional>
#include <vector>

namespace fmindex_collection {

struct CSA_SAS {
    std::vector<uint64_t> ssa;
    uint64_t samplingRate;

    CSA_SAS(std::vector<uint64_t> _ssa, uint64_t _samplingRate)
        : ssa{std::move(_ssa)}
        , samplingRate {_samplingRate}
    {}
    CSA_SAS(CSA_SAS const&) = delete;
    CSA_SAS(CSA_SAS&& _other) noexcept = default;
    CSA_SAS(cereal_tag) {}

    auto operator=(CSA_SAS const&) -> CSA_SAS& = delete;
    auto operator=(CSA_SAS&& _other) noexcept -> CSA_SAS& = default;

    size_t memoryUsage() const {
        return sizeof(ssa) + ssa.size() * sizeof(ssa.back());
    }

    auto value(size_t idx) const -> std::optional<uint64_t> {
        if (idx % samplingRate > 0) {
            return std::nullopt;
        }
        return ssa[idx / samplingRate];
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(ssa, samplingRate);
    }
};

}
