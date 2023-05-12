// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "concepts.h"

#include <array>
#include <cstdint>
#include <vector>

namespace fmindex_collection {
namespace occtable {
namespace naive {

template <uint64_t TSigma>
struct OccTable {
    using TLengthType = uint64_t;
    static constexpr uint64_t Sigma = TSigma;

    std::vector<std::array<uint64_t, Sigma>> occ{};
    std::array<uint64_t, Sigma+1> C{};

    static uint64_t expectedMemoryUsage(uint64_t length) {
        uint64_t C       = sizeof(uint64_t) * (Sigma+1);
        uint64_t entries = sizeof(uint64_t) * Sigma * (length+1);
        return C + entries;
    }

    OccTable(std::span<uint8_t const> _bwt) {
        occ.reserve(_bwt.size()+1);
        occ.push_back(std::array<uint64_t, Sigma>{});

        for (uint64_t i{0}; i < _bwt.size(); ++i) {
            occ.push_back(occ.back());
            occ.back()[_bwt[i]] += 1;
        }
        C[0] = 0;
        for (uint64_t i{0}; i < Sigma; ++i) {
            C[i+1] = occ.back()[i] + C[i];
        }
    }

    OccTable(cereal_tag) {}

    static auto name() -> std::string {
        return "Naive";
    }

    static auto extension() -> std::string {
        return "n";
    }

    uint64_t memoryUsage() const {
        return occ.size() * sizeof(occ.back()) + sizeof(OccTable);
    }

    uint64_t size() const {
        return C.back();
    }

    uint64_t rank(uint64_t idx, uint8_t symb) const {
        return occ[idx][symb] + C[symb];
    }

    uint64_t prefix_rank(uint64_t idx, uint8_t symb) const {
        uint64_t a{};
        for (uint64_t i{0}; i <= symb; ++i) {
            a += occ[idx][i];
        }
        return a;
    }

    uint64_t symbol(uint64_t idx) const {
        idx += 1;
        for (uint64_t i{0}; i < Sigma-1; ++i) {
            if (occ[idx][i] > occ[idx-1][i]) {
                return i;
            }
        }
        return Sigma-1;
    }

    auto all_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, Sigma>, std::array<uint64_t, Sigma>> {
        std::array<uint64_t, Sigma> rs{0};
        std::array<uint64_t, Sigma> prs{0};
        uint64_t a{0};
        for (uint64_t i{0}; i < Sigma; ++i) {
            auto r = occ[idx][i];
            rs[i] = r + C[i];
            a += r;
            prs[i] = a;
        }
        return {rs, prs};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(occ, C);
    }
};
static_assert(checkOccTable<OccTable>);

}
}
}
