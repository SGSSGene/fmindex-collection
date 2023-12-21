// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../builtins.h"
#include "../rankvector/rankvector.h"
#include "concepts.h"

#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <span>
#include <vector>

/**
 * A run length encoded bwt, see unpublished paper, but with recursive bitvectors
 *
 */
namespace fmindex_collection {
namespace occtable {
namespace generic_rrlebwt {


template <size_t Depth, size_t Sigma>
constexpr auto getNestedType() {
    if constexpr (Depth > 1) {
        return rankvector::RLE<Sigma, /*blockencodingsize*/2, rankvector::Naive<Sigma>, decltype(getNestedType<Depth-1, Sigma>())>{cereal_tag{}};
    } else {
        return rankvector::Naive<Sigma>{cereal_tag{}};
    }
}

template <uint64_t TSigma, uint64_t TAlignment, typename block_t>
struct OccTable {
    using TLengthType = uint64_t;
    static constexpr uint64_t Sigma = TSigma;

    using RankVector = decltype(getNestedType</*depth*/3, Sigma>());

    RankVector bitvector;
    std::array<size_t, TSigma+1> C{};

    static uint64_t expectedMemoryUsage(uint64_t length) {
        return 0;
    }

    OccTable(std::span<uint8_t const> _bwt)
        : bitvector{_bwt}
    {
        for (auto c : _bwt) {
            C[c+1] += 1;
        }
        for (size_t i{1}; i < C.size(); ++i) {
            C[i] = C[i] + C[i-1];
        }
    }

    OccTable(cereal_tag)
        : bitvector{cereal_tag{}}
    {}

    uint64_t memoryUsage() const {
        return 0;
    }

    auto prefetch(uint64_t idx) const {
        bitvector.prefetch(idx);
    }

    uint64_t size() const {
        return bitvector.size();
    }

    uint64_t symbol(uint64_t idx) const {
        return bitvector.symbol(idx);
    }

    uint64_t rank(uint64_t idx, uint64_t symb) const {
        return bitvector.rank(idx, symb) + C[symb];
    }

    uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
        return bitvector.prefix_rank(idx, symb);
    }


    auto all_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, Sigma>, std::array<uint64_t, Sigma>> {
        auto [rs, prs] = bitvector.all_ranks_and_prefix_ranks(idx);
        for (size_t i{0}; i < rs.size(); ++i) {
            rs[i] += C[i];
        }
        return {rs, prs};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bitvector, C);
    }
};
}

namespace rrlebwt {

template <uint64_t TSigma>
struct OccTable : generic_rrlebwt::OccTable<TSigma, 8, uint16_t> {
    using generic_rrlebwt::OccTable<TSigma, 8, uint16_t>::OccTable;
    static auto name() -> std::string {
        return "Interleaved 16bit, rrle";
    }

    static auto extension() -> std::string {
        return "i16rrle";
    }
};
static_assert(checkOccTable<OccTable>);
}

}
}
